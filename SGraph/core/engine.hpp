#ifndef ENGINE_HPP
#define ENGINE_HPP

#include <stdio.h>
#include <stdlib.h>

#include "core/graph.hpp"

enum ThreadStatus {
	WORKING,
	STEALING
};

struct ThreadState {
	VertexId curr;
	VertexId end;
	ThreadStatus status;
};

struct MessageBuffer {
	size_t capacity;
	int count; // the actual size (i.e. bytes) should be sizeof(element) * count
	char * data;
	MessageBuffer () {
		capacity = 0;
		count = 0;
		data = NULL;
	}
	void init () {
		capacity = 4096;
		count = 0;
		data = (char*)malloc(capacity);
	}
	void resize(size_t new_capacity) {
		if (new_capacity > capacity) {
			char * new_data = (char*)realloc(data, new_capacity);
			assert(new_data!=NULL);
			data = new_data;
			capacity = new_capacity;
		}
	}
};

template <typename MsgData>
struct MsgUnit {
	VertexId vertex;
	MsgData msg_data;
} __attribute__((packed));

template<typename EdgeData>
class Engine{
public:
	using Adjlist = std::vector<AdjEdge<EdgeData>>;

	Graph<EdgeData> * graph;
    size_t alpha;
	int threads;
    int message_tag;
	MPI_Comm& comm;	
    
    int partition_id;
	int partitions;

	VertexId * partition_offset; // VertexId [partitions+1]
	VertexId owned_vertices;

	VertexId * out_degree;	// VertexId [vertices];
	VertexId * in_degree; // VertexId [vertices];	

	Storage<EdgeData> ** outgoing_storage; // Storage * [vertices]; 
	Storage<EdgeData> ** incoming_storage; // Storage * [vertices]; 

	ThreadState ** thread_state; // ThreadState* [threads]; 
	ThreadState ** tuned_chunks_dense; // ThreadState [partitions][threads];
	ThreadState ** tuned_chunks_sparse; // ThreadState [partitions][threads];

	size_t local_send_buffer_limit;
	MessageBuffer ** local_send_buffer; // MessageBuffer* [threads];

	int current_send_part_id;
	MessageBuffer ** send_buffer; // MessageBuffer* [partitions];
	MessageBuffer ** recv_buffer; // MessageBuffer* [partitions];

    Engine(Graph<EdgeData>* _graph, int _threads, int _message_tag, MPI_Comm& _comm) : comm(_comm) {
        graph = _graph;
        partition_id = _graph->partition_id;
        partitions = _graph->partitions;
        partition_offset = graph->partition_offset;
        owned_vertices = graph->owned_vertices;
        out_degree = graph->out_degree;
        in_degree = graph->in_degree;
        outgoing_storage = graph->outgoing_storage;
        incoming_storage = graph->incoming_storage;
	    alpha = 8 * (graph->partitions - 1) + 1;
        threads = _threads;
        message_tag = _message_tag;

		thread_state = new ThreadState * [threads];
		local_send_buffer_limit = 16;
		local_send_buffer = new MessageBuffer * [threads];
		for (int t_i=0;t_i<threads;t_i++) {
			thread_state[t_i] = (ThreadState*)malloc( sizeof(ThreadState));
			local_send_buffer[t_i] = (MessageBuffer*)malloc( sizeof(MessageBuffer));
			local_send_buffer[t_i]->init();
		}

		send_buffer = new MessageBuffer * [partitions];
		recv_buffer = new MessageBuffer * [partitions];
		for (int i=0;i<partitions;i++) {
			send_buffer[i] = (MessageBuffer*)malloc( sizeof(MessageBuffer));
			recv_buffer[i] = (MessageBuffer*)malloc( sizeof(MessageBuffer));
			send_buffer[i]->init();
			recv_buffer[i]->init();
		}

		transpose();
		tune_chunks();
		transpose();
		tune_chunks();
    }

	// transpose the graph
	void transpose() {
		std::swap(out_degree, in_degree);
		std::swap(outgoing_storage, incoming_storage);
		std::swap(tuned_chunks_dense, tuned_chunks_sparse);
	}

	void tune_chunks() {
		size_t basic_chunk = 64;
		tuned_chunks_dense = new ThreadState * [partitions];
//        int iter = 0;
		for (int p_i=0;p_i<partitions;p_i++) {
			tuned_chunks_dense[p_i] = new ThreadState [threads];
			int remained_threads = threads;
			EdgeId remained_size = 0;
			for (VertexId v_i=partition_offset[p_i];v_i<partition_offset[p_i+1];v_i++) {
				// remained_size += graph->incoming_storage[v_i]->adjlist_size;
				remained_size += 1;
				remained_size += alpha;
//                iter++;
			}
//            printf("remained_size: %d\n", remained_size);
//            printf("iter: %d\n", iter);
			VertexId last_v_i = partition_offset[p_i];
			for (int t_i=0;t_i<threads;t_i++){
				tuned_chunks_dense[p_i][t_i].status = WORKING;
				tuned_chunks_dense[p_i][t_i].curr = last_v_i;
				if (remained_threads == 1) {
					tuned_chunks_dense[p_i][t_i].end = partition_offset[p_i+1];
				}
				else {
					EdgeId expected_size = remained_size / remained_threads;
//                    printf("expected_size: %d\n", expected_size);
					tuned_chunks_dense[p_i][t_i].end = last_v_i;
					for (VertexId v_i=last_v_i;v_i+basic_chunk<partition_offset[p_i+1];v_i+=basic_chunk) {
						EdgeId got_size = 0;
						for (VertexId v_j=v_i;v_j<v_i+basic_chunk;v_j++) {
							// got_size += graph->incoming_storage[v_j]->adjlist_size;
							got_size += 1;
							got_size += alpha;
						}
						if (got_size <= expected_size) {
							tuned_chunks_dense[p_i][t_i].end += basic_chunk;
							expected_size -= got_size;
							remained_size -= got_size;
						}
						else {
							break;
						}
					}
					last_v_i = tuned_chunks_dense[p_i][t_i].end;
					remained_threads -= 1;
				}
			}
		}
	}	
    
    // process vertices
	template<typename R>
	R process_vertices(std::function<R(VertexId)> process, Bitmap * active, bool trace=false) {
		double stream_time = 0;
		stream_time -= MPI_Wtime();

		R reducer = 0;
		size_t basic_chunk = 64;
		for (int t_i=0;t_i<threads;t_i++) {
			thread_state[t_i]->curr = partition_offset[partition_id] + owned_vertices / threads * t_i / basic_chunk * basic_chunk;
			thread_state[t_i]->end = partition_offset[partition_id] + owned_vertices / threads * (t_i+1)  / basic_chunk * basic_chunk;
			if (t_i == threads - 1) {
				thread_state[t_i]->end = partition_offset[partition_id+1];
			}
			thread_state[t_i]->status = WORKING;
		}
		#pragma omp parallel reduction(+:reducer)
		{
			R local_reducer = 0;
			int thread_id = omp_get_thread_num();
			while (true) {
				VertexId v_i = __sync_fetch_and_add(&thread_state[thread_id]->curr, basic_chunk);
				if (v_i >= thread_state[thread_id]->end) break;
				unsigned long word = active->data[WORD_OFFSET(v_i)];
				while (word != 0) {
					if (word & 1) {
                        if(trace == true) printf("Thread [%d] process vertex [%ld]\n", thread_id, v_i);
						local_reducer += process(v_i);
					}
					v_i++;
					word = word >> 1;
				}
			}
			thread_state[thread_id]->status = STEALING;
			for (int t_offset=1;t_offset<threads;t_offset++) {
				int t_i = (thread_id + t_offset) % threads;
				while (thread_state[t_i]->status!=STEALING) {
					VertexId v_i = __sync_fetch_and_add(&thread_state[t_i]->curr, basic_chunk);
					if (v_i >= thread_state[t_i]->end) continue;
					unsigned long word = active->data[WORD_OFFSET(v_i)];
					while (word != 0) {
						if (word & 1) {
                            if(trace == true) printf("Thread [%d] process vertex [%ld]\n", thread_id, v_i);
							local_reducer += process(v_i);
						}
						v_i++;
						word = word >> 1;
					}
				}
			}
			reducer += local_reducer;
		}
		R global_reducer;
		MPI_Datatype dt = get_mpi_data_type<R>();
		MPI_Allreduce(&reducer, &global_reducer, 1, dt, MPI_SUM, comm);
		stream_time += MPI_Wtime();
		#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id==0) {
			printf("process_vertices took %lf (s)\n", stream_time);
		}
		#endif
		return global_reducer;
	}

	template<typename M>
	void flush_local_send_buffer(int t_i) {
		int pos = __sync_fetch_and_add(&send_buffer[current_send_part_id]->count, local_send_buffer[t_i]->count);
		memcpy(send_buffer[current_send_part_id]->data + sizeof(MsgUnit<M>) * pos, local_send_buffer[t_i]->data, sizeof(MsgUnit<M>) * local_send_buffer[t_i]->count);
		local_send_buffer[t_i]->count = 0;
	}

	// emit a message to a vertex's master (dense) / mirror (sparse)
	template<typename M>
	void emit(VertexId vtx, M msg) {
    	int t_i = omp_get_thread_num();
		MsgUnit<M> * buffer = (MsgUnit<M>*)local_send_buffer[t_i]->data;
		buffer[local_send_buffer[t_i]->count].vertex = vtx;
		buffer[local_send_buffer[t_i]->count].msg_data = msg;
		local_send_buffer[t_i]->count += 1;
		if (local_send_buffer[t_i]->count==local_send_buffer_limit) {
			flush_local_send_buffer<M>(t_i);
		}
	}

	template<typename R, typename M>
	R process_edges_sparse(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, Bitmap * active) {
		for (int t_i=0;t_i<threads;t_i++) {
			local_send_buffer[t_i]->resize( sizeof(MsgUnit<M>) * local_send_buffer_limit );
			local_send_buffer[t_i]->count = 0;
		}
		double stream_time = 0;
		stream_time -= MPI_Wtime();
		R reducer = 0;
		for (int i=0;i<partitions;i++) {
			recv_buffer[i]->resize( sizeof(MsgUnit<M>) * (partition_offset[i+1] - partition_offset[i]) );
			send_buffer[i]->resize( sizeof(MsgUnit<M>) * owned_vertices );
			send_buffer[i]->count = 0;
			recv_buffer[i]->count = 0;
		}
		size_t basic_chunk = 64;
		#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id==0) {
			printf("sparse mode\n");
		}
		#endif
		int * recv_queue = new int [partitions];
		int recv_queue_size = 0;
		std::mutex recv_queue_mutex;

		current_send_part_id = partition_id;
		#pragma omp parallel for
		for (VertexId begin_v_i=partition_offset[partition_id];begin_v_i<partition_offset[partition_id+1];begin_v_i+=basic_chunk) {
			VertexId v_i = begin_v_i;
			unsigned long word = active->data[WORD_OFFSET(v_i)];
			while (word != 0) {
				if (word & 1) {
					sparse_signal(v_i);
				}
				v_i++;
				word = word >> 1;
			}
		}

		#pragma omp parallel for
		for (int t_i=0;t_i<threads;t_i++) {
			flush_local_send_buffer<M>(t_i);
		}
		recv_queue[recv_queue_size] = partition_id;
		recv_queue_mutex.lock();
		recv_queue_size += 1;
		recv_queue_mutex.unlock();
		std::thread send_thread([&](){
			for (int step=1;step<partitions;step++) {
				int i = (partition_id - step + partitions) % partitions;
				MPI_Send(send_buffer[partition_id]->data, sizeof(MsgUnit<M>) * send_buffer[partition_id]->count, MPI_CHAR, i, message_tag, comm);
			}
		});
		std::thread recv_thread([&](){
			for (int step=1;step<partitions;step++) {
				int i = (partition_id + step) % partitions;
				MPI_Status recv_status;
				MPI_Probe(i, message_tag, comm, &recv_status);
				MPI_Get_count(&recv_status, MPI_CHAR, &recv_buffer[i]->count);
				MPI_Recv(recv_buffer[i]->data, recv_buffer[i]->count, MPI_CHAR, i, message_tag, comm, MPI_STATUS_IGNORE);
				recv_buffer[i]->count /= sizeof(MsgUnit<M>);
				recv_queue[recv_queue_size] = i;
				recv_queue_mutex.lock();
				recv_queue_size += 1;
				recv_queue_mutex.unlock();
			}
		});
		for (int step=0;step<partitions;step++) {
			while (true) {
				recv_queue_mutex.lock();
				bool condition = (recv_queue_size<=step);
				recv_queue_mutex.unlock();
				if (!condition) break;
				__asm volatile ("pause" ::: "memory");
			}
			int i = recv_queue[step];
			MessageBuffer * used_buffer;
			if (i==partition_id) {
				used_buffer = send_buffer[i];
			} else {
				used_buffer = recv_buffer[i];
			}
			MsgUnit<M> * buffer = (MsgUnit<M> *)used_buffer->data;
			size_t buffer_size = used_buffer->count;

			// test omp vs work stealing
			// #pragma omp parallel for reduction(+:reducer)
			// for (size_t b_i=0;b_i<buffer_size;b_i++)
			// {
			// 	VertexId v_i = buffer[b_i].vertex;
			// 	M msg_data = buffer[b_i].msg_data;
			// 	reducer += sparse_slot(v_i, msg_data, outgoing_adjlist[v_i]);	
			// }

			for (int t_i=0;t_i<threads;t_i++) {
				VertexId partition_size = buffer_size;
				thread_state[t_i]->curr = partition_size / threads * t_i / basic_chunk * basic_chunk;
				thread_state[t_i]->end = partition_size /threads * (t_i+1) / basic_chunk * basic_chunk;
				if (t_i == threads - 1) {
					thread_state[t_i]->end = buffer_size;
				}
				thread_state[t_i]->status = WORKING;
			}
			#pragma omp parallel reduction(+:reducer)
			{
				R local_reducer = 0;
				int thread_id = omp_get_thread_num();
				while (true) {
					VertexId b_i = __sync_fetch_and_add(&thread_state[thread_id]->curr, basic_chunk);
					if (b_i >= thread_state[thread_id]->end) break;
					VertexId begin_b_i = b_i;
					VertexId end_b_i = b_i + basic_chunk;
					if (end_b_i>thread_state[thread_id]->end) {
						end_b_i = thread_state[thread_id]->end;
					}
					for (b_i=begin_b_i;b_i<end_b_i;b_i++) {
						VertexId v_i = buffer[b_i].vertex;
						M msg_data = buffer[b_i].msg_data;
						local_reducer += sparse_slot(v_i, msg_data, outgoing_storage[v_i]->adjlist);
					}
				}
				thread_state[thread_id]->status = STEALING;
				for (int t_offset=1;t_offset<threads;t_offset++) {
					int t_i = (thread_id + t_offset) % threads;
					if (thread_state[t_i]->status==STEALING) continue;
					while (true) {
						VertexId b_i = __sync_fetch_and_add(&thread_state[t_i]->curr, basic_chunk);
						if (b_i >= thread_state[t_i]->end) break;
						VertexId begin_b_i = b_i;
						VertexId end_b_i = b_i + basic_chunk;
						if (end_b_i>thread_state[t_i]->end) {
							end_b_i = thread_state[t_i]->end;
						}
						for (b_i=begin_b_i;b_i<end_b_i;b_i++) {
							VertexId v_i = buffer[b_i].vertex;
							M msg_data = buffer[b_i].msg_data;
							local_reducer += sparse_slot(v_i, msg_data, outgoing_storage[v_i]->adjlist);
						}
					}
				}
				reducer += local_reducer;
			}
		}
		send_thread.join();
		recv_thread.join();
		delete [] recv_queue;

		R global_reducer;
		MPI_Datatype dt = get_mpi_data_type<R>();
		MPI_Allreduce(&reducer, &global_reducer, 1, dt, MPI_SUM, comm);
		stream_time += MPI_Wtime();
		#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id==0) {
			printf("process_edges took %lf (s)\n", stream_time);
		}
		#endif
		return global_reducer;
	}

	template<typename R, typename M>
	R process_edges_dense(std::function<void(VertexId, Adjlist&)> dense_signal, std::function<R(VertexId, M)> dense_slot, Bitmap * active, Bitmap * dense_selective = nullptr, Bitmap ** dense_selective_hubs = nullptr) {
		for (int t_i=0;t_i<threads;t_i++) {
			local_send_buffer[t_i]->resize( sizeof(MsgUnit<M>) * local_send_buffer_limit );
			local_send_buffer[t_i]->count = 0;
		}
		double stream_time = 0;
		stream_time -= MPI_Wtime();
		R reducer = 0;
		for (int i=0;i<partitions;i++) {
			recv_buffer[i]->resize( sizeof(MsgUnit<M>) * owned_vertices );
			send_buffer[i]->resize( sizeof(MsgUnit<M>) * (partition_offset[i+1] - partition_offset[i]) );
			send_buffer[i]->count = 0;
			recv_buffer[i]->count = 0;
		}
		size_t basic_chunk = 64;
		// dense selective bitmap
		if (dense_selective!=nullptr && partitions>1) {
			double sync_time = 0;
			sync_time -= get_time();
			std::thread send_thread([&](){
				for (int step=1;step<partitions;step++) {
					int recipient_id = (partition_id + step) % partitions;
					VertexId words = (partition_offset[partition_id+1] - partition_offset[partition_id]) / 64;
					if ((partition_offset[partition_id+1] - partition_offset[partition_id]) % 64) {
						words++;
					}
					MPI_Send(dense_selective->data + WORD_OFFSET(partition_offset[partition_id]), words, MPI_UNSIGNED_LONG, recipient_id, message_tag, comm);
					if (dense_selective_hubs != nullptr) {
						for (int i=0;i<hubs;i++) {
							MPI_Send(dense_selective_hubs[i]->data + WORD_OFFSET(partition_offset[partition_id]), words, MPI_UNSIGNED_LONG, recipient_id, message_tag, comm);
						}
					}
				}
			});
			std::thread recv_thread([&](){
				for (int step=1;step<partitions;step++) {
					int sender_id = (partition_id - step + partitions) % partitions;
					VertexId words = (partition_offset[sender_id+1] - partition_offset[sender_id]) / 64;
					if ((partition_offset[sender_id+1] - partition_offset[sender_id]) % 64) {
						words++;
					}
					MPI_Recv(dense_selective->data + WORD_OFFSET(partition_offset[sender_id]), words, MPI_UNSIGNED_LONG, sender_id, message_tag, comm, MPI_STATUS_IGNORE);
					if (dense_selective_hubs != nullptr) {
						for (int i=0;i<hubs;i++) {
							MPI_Recv(dense_selective_hubs[i]->data + WORD_OFFSET(partition_offset[sender_id]), words, MPI_UNSIGNED_LONG, sender_id, message_tag, comm, MPI_STATUS_IGNORE);
						}
					}
				}
			});
			send_thread.join();
			recv_thread.join();
			sync_time += get_time();
			#ifdef PRINT_DEBUG_MESSAGES
			if (partition_id==0) {
				printf("sync_time = %lf\n", sync_time);
			}
			#endif
		}
		#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id==0) {
			printf("dense mode\n");
		}
		#endif
		int * send_queue = new int [partitions];
		int * recv_queue = new int [partitions];
		volatile int send_queue_size = 0;
		volatile int recv_queue_size = 0;
		std::mutex send_queue_mutex;
		std::mutex recv_queue_mutex;

		std::thread send_thread([&](){
			for (int step=0;step<partitions;step++) {
				if (step==partitions-1) {
					break;
				}
				while (true) {
					send_queue_mutex.lock();
					bool condition = (send_queue_size<=step);
					send_queue_mutex.unlock();
					if (!condition) break;
					__asm volatile ("pause" ::: "memory");
				}
				int i = send_queue[step];
				MPI_Send(send_buffer[i]->data, sizeof(MsgUnit<M>) * send_buffer[i]->count, MPI_CHAR, i, message_tag, comm);
			}
		});
		std::thread recv_thread([&](){
			std::vector<std::thread> threads;
			for (int step=1;step<partitions;step++) {
				int i = (partition_id - step + partitions) % partitions;
				threads.emplace_back([&](int i){
					MPI_Status recv_status;
					MPI_Probe(i, message_tag, comm, &recv_status);
					MPI_Get_count(&recv_status, MPI_CHAR, &recv_buffer[i]->count);
					MPI_Recv(recv_buffer[i]->data, recv_buffer[i]->count, MPI_CHAR, i, message_tag, comm, MPI_STATUS_IGNORE);
					recv_buffer[i]->count /= sizeof(MsgUnit<M>);
				}, i);
			}
			for (int step=1;step<partitions;step++) {
				int i = (partition_id - step + partitions) % partitions;
				threads[step-1].join();
				recv_queue[recv_queue_size] = i;
				recv_queue_mutex.lock();
				recv_queue_size += 1;
				recv_queue_mutex.unlock();
			}
			recv_queue[recv_queue_size] = partition_id;
			recv_queue_mutex.lock();
			recv_queue_size += 1;
			recv_queue_mutex.unlock();
		});
		current_send_part_id = partition_id;
		for (int step=0;step<partitions;step++) {
			current_send_part_id = (current_send_part_id + 1) % partitions;
			int i = current_send_part_id;
			for (int t_i=0;t_i<threads;t_i++) {
				*thread_state[t_i] = tuned_chunks_dense[i][t_i];
			}
			#pragma omp parallel
			{
				int thread_id = omp_get_thread_num();
				VertexId final_p_v_i = thread_state[thread_id]->end;
				while (true) {
					VertexId begin_p_v_i = __sync_fetch_and_add(&thread_state[thread_id]->curr, basic_chunk);
					if (begin_p_v_i >= final_p_v_i) break;
					VertexId end_p_v_i = begin_p_v_i + basic_chunk;
					if (end_p_v_i > final_p_v_i) {
						end_p_v_i = final_p_v_i;
					}
					for (VertexId p_v_i = begin_p_v_i; p_v_i < end_p_v_i; p_v_i ++) {
						VertexId v_i = p_v_i;
						dense_signal(v_i, incoming_storage[v_i]->adjlist);
					}
				}
				thread_state[thread_id]->status = STEALING;
				for (int t_offset=1;t_offset<threads;t_offset++) {
					int t_i = (thread_id + t_offset) % threads;
					while (thread_state[t_i]->status!=STEALING) {
						VertexId begin_p_v_i = __sync_fetch_and_add(&thread_state[t_i]->curr, basic_chunk);
						if (begin_p_v_i >= thread_state[t_i]->end) break;
						VertexId end_p_v_i = begin_p_v_i + basic_chunk;
						if (end_p_v_i > thread_state[t_i]->end) {
							end_p_v_i = thread_state[t_i]->end;
						}
						for (VertexId p_v_i = begin_p_v_i; p_v_i < end_p_v_i; p_v_i ++) {
							VertexId v_i = p_v_i;
							dense_signal(v_i, incoming_storage[v_i]->adjlist);
						}
					}
				}
			}
			#pragma omp parallel for
			for (int t_i=0;t_i<threads;t_i++) {
				flush_local_send_buffer<M>(t_i);
			}
			if (i!=partition_id) {
				send_queue[send_queue_size] = i;
				send_queue_mutex.lock();
				send_queue_size += 1;
				send_queue_mutex.unlock();
			}
		}
		for (int step=0;step<partitions;step++) {
			while (true) {
				recv_queue_mutex.lock();
				bool condition = (recv_queue_size<=step);
				recv_queue_mutex.unlock();
				if (!condition) break;
				__asm volatile ("pause" ::: "memory");
			}
			int i = recv_queue[step];
			MessageBuffer * used_buffer;
			if (i==partition_id) {
				used_buffer = send_buffer[i];
			} else {
				used_buffer = recv_buffer[i];
			}
			for (int t_i=0;t_i<threads;t_i++) {
				VertexId partition_size = used_buffer->count;
				thread_state[t_i]->curr = partition_size / threads * t_i / basic_chunk * basic_chunk;
				thread_state[t_i]->end = partition_size /threads * (t_i+1) / basic_chunk * basic_chunk;
				if (t_i == threads - 1) {
					thread_state[t_i]->end = used_buffer->count;
				}
				thread_state[t_i]->status = WORKING;
			}
			#pragma omp parallel reduction(+:reducer)
			{
				R local_reducer = 0;
				int thread_id = omp_get_thread_num();
				MsgUnit<M> * buffer = (MsgUnit<M> *)used_buffer->data;
				while (true) {
					VertexId b_i = __sync_fetch_and_add(&thread_state[thread_id]->curr, basic_chunk);
					if (b_i >= thread_state[thread_id]->end) break;
					VertexId begin_b_i = b_i;
					VertexId end_b_i = b_i + basic_chunk;
					if (end_b_i>thread_state[thread_id]->end) {
						end_b_i = thread_state[thread_id]->end;
					}
					for (b_i=begin_b_i;b_i<end_b_i;b_i++) {
						VertexId v_i = buffer[b_i].vertex;
						M msg_data = buffer[b_i].msg_data;
						local_reducer += dense_slot(v_i, msg_data);
					}
				}
				thread_state[thread_id]->status = STEALING;
				reducer += local_reducer;
			}
		}
		send_thread.join();
		recv_thread.join();
		delete [] send_queue;
		delete [] recv_queue;

		R global_reducer;
		MPI_Datatype dt = get_mpi_data_type<R>();
		MPI_Allreduce(&reducer, &global_reducer, 1, dt, MPI_SUM, comm);
		stream_time += MPI_Wtime();
		#ifdef PRINT_DEBUG_MESSAGES
		if (partition_id==0) {
			printf("process_edges took %lf (s)\n", stream_time);
		}
		#endif
		return global_reducer;
	}

	// process edges
	template<typename R, typename M>
	R process_edges(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, std::function<void(VertexId, Adjlist&)> dense_signal, std::function<R(VertexId, M)> dense_slot, Bitmap * active, Bitmap * dense_selective = nullptr, Bitmap ** dense_selective_hubs = nullptr) {
		EdgeId active_edges = process_vertices<EdgeId>(
			[&](VertexId vtx){
				return (EdgeId)out_degree[vtx];
			},
			active
		);
		bool sparse = (active_edges < graph->edge_num * 0.05);
		// bool sparse = true;
		if (sparse) return process_edges_sparse(sparse_signal, sparse_slot, active);
		else return process_edges_dense(dense_signal, dense_slot, active, dense_selective, dense_selective_hubs);
	}
};

#endif
