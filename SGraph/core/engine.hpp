#ifndef ENGINE_HPP
#define ENGINE_HPP

#include <stdio.h>
#include <stdlib.h>

#include "core/graph.hpp"

enum ThreadStatus {//定义了一个枚举类型ThreadStatus，用于表示线程的状态:WORKING表示线程正在工作，STEALING表示线程正在窃取其他线程的任务。
	WORKING,
	STEALING
};

struct ThreadState {//定义了一个结构体ThreadState，用于表示线程的状态。该结构体包含了当前处理的顶点的起始和结束索引，以及线程的状态。
	VertexId curr;
	VertexId end;
	ThreadStatus status;
};

struct MessageBuffer {//定义了一个结构体MessageBuffer，用于表示消息缓冲区。该结构体包含了缓冲区的容量、当前的大小以及缓冲区的数据。
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
	int threads;//每个分区的线程数
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
	    alpha = 8 * (graph->partitions - 1) + 1;//作者说，这里的设置其实没有多大的影响，反正每个顶点都是加上一样的值
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
		tune_chunks();//据作者说，这里可以省略
    }

	// transpose the graph
	void transpose() {
		std::swap(out_degree, in_degree);
		std::swap(outgoing_storage, incoming_storage);
		std::swap(tuned_chunks_dense, tuned_chunks_sparse);
	}


	//这个函数特别绕，因为里面牵涉到作者对工作负载的理解，以及根据工作负载为每个线程分配顶点数的逻辑。
	void tune_chunks() {//Engine 类中的多线程操作需要对图中的顶点进行分块处理，以便于每个线程能够操作一个特定的顶点块。tune_chunks 函数的目标就是使每个线程处理的负载尽可能均衡。
		size_t basic_chunk = 64;//定义了基本的块大小，即每个块包含的顶点数。
		tuned_chunks_dense = new ThreadState * [partitions];//tuned_chunks_dense: 一个二维数组，用于存储每个分区和线程的状态。


		////////////统计总工作量
		for (int p_i=0;p_i<partitions;p_i++) {//遍历每一个分区 p_i，对每个分区都进行调整
			tuned_chunks_dense[p_i] = new ThreadState [threads];//初始化 tuned_chunks_dense 为每个分区的线程状态的数组。
			int remained_threads = threads;
			EdgeId remained_size = 0;//对于当前分区，初始化 remained_threads 为总的线程数和 remained_size 为0。
			for (VertexId v_i=partition_offset[p_i];v_i<partition_offset[p_i+1];v_i++) {//这个循环遍历分区内所有的顶点，将每个顶点的工作量都累加，计算得到总的工作量remained_size
				// remained_size += graph->incoming_storage[v_i]->adjlist_size;
				//这里注释了一行代码，原始的gemini在计算每个顶点的工作量时，把每个顶点的度数累加得到各自的工作量。但是作者对原始图数据集进行了一次shuffle，所以只要分区的顶点数一样，就认为工作量大致相当。
				remained_size += 1;
				remained_size += alpha;// 为当前分区中的每个顶点计算 remained_size。这里，remained_size 对于每个顶点增加1和 alpha 值。
				/*
					remained_size代表待处理的数据量。 这个变量名"remained"意味着它可能表示待处理或尚未处理的部分。因为它在循环中逐渐增加，它很可能表示某种累积的数据量。
					remained_size与待处理的顶点的数量有关。 这从循环可以看出，它遍历了一个分区中的所有顶点，累计所有顶点的空间
					每个顶点的度数都不同，理论上remained_size等于待处理的顶点的度数之和。现在为了方便计算，直接默认每个顶点的度数为alpha，所以remained_size每次加1加alpha。 
				*/
			}


			////////////////根据上面得到的工作量，尽可能为每个线程均分工作负载
			VertexId last_v_i = partition_offset[p_i];//使用一个 last_v_i 变量来追踪当前分区的最后一个处理的顶点。
			for (int t_i=0;t_i<threads;t_i++){//为当前分区的每个线程 t_i 分配一个块：

				//准备处理线程t_i
				tuned_chunks_dense[p_i][t_i].status = WORKING;//设置该线程的状态为 `WORKING`。
				tuned_chunks_dense[p_i][t_i].curr = last_v_i;//为分区 p_i 中的线程 t_i 设置当前处理的块的起始顶点为 last_v_i。在后续的循环中，last_v_i 会被更新为下一个块的起始顶点。
				
				//如果剩余1个线程（线程t_i是最后一个线程），剩下的所有顶点都分给该线程
				if (remained_threads == 1) {//设定块的起始和结束顶点。这里的目标是尽量使每个线程处理的顶点数量均等。但是，如果当前是最后一个线程，则确保处理分区内的所有剩余顶点。
					tuned_chunks_dense[p_i][t_i].end = partition_offset[p_i+1];
				}
				else {//如果还有多于一个线程，计算期望的大小（`expected_size`），这是为了确保每个线程处理大致相同的工作量
					EdgeId expected_size = remained_size / remained_threads;
					tuned_chunks_dense[p_i][t_i].end = last_v_i;
					for (VertexId v_i=last_v_i;v_i+basic_chunk<partition_offset[p_i+1];v_i+=basic_chunk) {//在一个分区内，以basic_chunk为单位处理顶点，直到达到分区的末尾或接近末尾。
						//一个顶点对应的工作量是（1+alpha），一次处理一个块，对应的工作量是（basic_chunk*（1+alpha））。要求工作量尽可能大，但是要小于expected_size
						EdgeId got_size = 0;
						for (VertexId v_j=v_i;v_j<v_i+basic_chunk;v_j++) {
							// got_size += graph->incoming_storage[v_j]->adjlist_size;
							got_size += 1;
							got_size += alpha;//这里其实很没有意义，因为作者在对图数据集进行shuffle后，所有的顶点被一视同仁，每个顶点的工作量都是（1+alpha），不如把每个顶点的工作量设为1。但是我们没有对图数据集进行shuffle，所以可以累加顶点的度数来计算工作负载。
						}
						if (got_size <= expected_size) {
							tuned_chunks_dense[p_i][t_i].end += basic_chunk;
							expected_size -= got_size;
							remained_size -= got_size;
						}
						else {//如果工作量已经大于expected_size，这时不再为线程分配新的顶点。
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
	R process_vertices(std::function<R(VertexId)> process, Bitmap * active) {//该函数的目的是并行地处理活动的顶点，并累积每个顶点处理的结果。处理完成后，所有线程的结果会被合并，并返回这个聚合的结果。
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
						local_reducer += process(v_i);
					}
					v_i++;
					word = word >> 1;
				}
			}
			thread_state[thread_id]->status = STEALING;
			for (int t_offset=1;t_offset<threads;t_offset++) {//当线程完成其分配的任务后，它会尝试窃取其他线程的任务，即处理其他线程尚未处理的顶点。此“窃取”策略有助于确保所有线程均衡地进行工作，并能够提高效率。
				int t_i = (thread_id + t_offset) % threads;
				while (thread_state[t_i]->status!=STEALING) {//这个循环会处理分配给当前线程的所有顶点。通过__sync_fetch_and_add函数来原子性地更新curr并获取当前顶点索引v_i。如果这个索引超过了该线程应该处理的范围，循环就会终止。
					VertexId v_i = __sync_fetch_and_add(&thread_state[t_i]->curr, basic_chunk);
					if (v_i >= thread_state[t_i]->end) continue;
					unsigned long word = active->data[WORD_OFFSET(v_i)];
					while (word != 0) {
						if (word & 1) {
							local_reducer += process(v_i);//process函数会使用“上界剪枝”+“下界剪枝”来判断顶点是否需要处理，如果需要处理就返回1，如果会被剪枝就返回0。所以这里的reducer存放的应该是需要处理的顶点的数量。
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
	void flush_local_send_buffer(int t_i) {//将本地的发送缓冲区的数据复制到全局的发送缓冲区，并重置本地发送缓冲区的计数。
		int pos = __sync_fetch_and_add(&send_buffer[current_send_part_id]->count, local_send_buffer[t_i]->count);
		memcpy(send_buffer[current_send_part_id]->data + sizeof(MsgUnit<M>) * pos, local_send_buffer[t_i]->data, sizeof(MsgUnit<M>) * local_send_buffer[t_i]->count);
		//使用memcpy函数将local_send_buffer[t_i]中的数据复制到send_buffer[current_send_part_id]的相应位置。	这里通过计算目标地址的偏移量（使用sizeof(MsgUnit<M>) * pos）来确定数据应该复制到send_buffer的哪个位置。
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
			flush_local_send_buffer<M>(t_i);//如果当前线程的本地发送缓冲区已满（即达到了local_send_buffer_limit的限制），则调用flush_local_send_buffer<M>(t_i)函数来清空缓冲区。这是为了确保缓冲区不会溢出，并且消息可以及时发送。	
		}
	}

	template<typename R, typename M>
	R process_edges_sparse(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, Bitmap * active) {
		for (int t_i=0;t_i<threads;t_i++) {//对每个线程，都设置了一个局部的消息发送缓冲区大小并将计数器设为0
			local_send_buffer[t_i]->resize( sizeof(MsgUnit<M>) * local_send_buffer_limit );
			local_send_buffer[t_i]->count = 0;
		}
		double stream_time = 0;
		stream_time -= MPI_Wtime();
		R reducer = 0;
		for (int i=0;i<partitions;i++) {//和前一个缓冲区不同，这里是分区对应的接收与发送缓冲区的初始化
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
			unsigned long word = active->data[WORD_OFFSET(v_i)];//WORD_OFFSET(v_i)是一个宏，用于确定给定顶点v_i在位图中的位置
			while (word != 0) {//这是一个while循环，用于遍历word中的所有位，直到word为0（即没有更多的活跃顶点）。
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
		recv_queue[recv_queue_size] = partition_id;//此时recv_queue_size为0，第一条消息的分区编号为partition_id，即当前分区的编号。
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
				int i = (partition_id + step) % partitions;//这一步很巧妙：自己当前所在的分区是partition_id，这一句通过step来遍历剩下的（partitions-1）个分区。同时用%partitions来保证得到的分区编号是有效的。
				MPI_Status recv_status;
				MPI_Probe(i, message_tag, comm, &recv_status);//MPI_Probe函数用于检查是否有消息到达。如果有消息到达，它会返回消息的大小和来源进程的编号。它的参数：来源进程的编号、消息的标签、通信域、消息状态。
				MPI_Get_count(&recv_status, MPI_CHAR, &recv_buffer[i]->count);//MPI_Get_count函数用于获取消息的大小。它的参数：消息状态、数据类型、消息大小。
				MPI_Recv(recv_buffer[i]->data, recv_buffer[i]->count, MPI_CHAR, i, message_tag, comm, MPI_STATUS_IGNORE);//MPI_Recv函数用于接收消息。它的参数：接收缓冲区、消息大小、数据类型、来源进程的编号、消息标签、通信域、消息状态。
				recv_buffer[i]->count /= sizeof(MsgUnit<M>);//这里将接收缓冲区的大小转换为消息单元的数量。
				recv_queue[recv_queue_size] = i;//recv_queue_size从0开始自增，它记录接收队列收到的消息分别来自哪个分区
				recv_queue_mutex.lock();
				recv_queue_size += 1;
				recv_queue_mutex.unlock();
			}
		});
		for (int step=0;step<partitions;step++) {
			while (true) {//确保recv_queue_size>step)
				recv_queue_mutex.lock();
				bool condition = (recv_queue_size<=step);
				recv_queue_mutex.unlock();
				if (!condition) break;
				__asm volatile ("pause" ::: "memory");//这里使用了一个内联汇编语句，它会告诉CPU暂停执行当前线程，直到有新的消息到达。
			}
			int i = recv_queue[step];
			MessageBuffer * used_buffer;
			if (i==partition_id) {//如果消息来自当前分区，则使用分区对应的发送缓冲区。
				used_buffer = send_buffer[i];
			} else {//如果消息来自其他分区，则使用分区对对应的接收缓冲区。
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

			for (int t_i=0;t_i<threads;t_i++) {//给当前分区的每个线程均分工作负载
				VertexId partition_size = buffer_size;//当前分区要处理的总工作负载
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
				while (true) {//每个线程分到自己的工作负载后，会以basic_chunk为单位处理，直到将所有负载处理完
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
						local_reducer += sparse_slot(v_i, msg_data, outgoing_storage[v_i]->adjlist);//sparse_slot函数会遍历出边更新上界，返回被激活的顶点数量。
					}
				}
				thread_state[thread_id]->status = STEALING;
				for (int t_offset=1;t_offset<threads;t_offset++) {
					int t_i = (thread_id + t_offset) % threads;//和这一句类似“int i = (partition_id + step) % partitions;” 主要目的是遍历所有除自己外的线程
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
		send_thread.join();//等待发送线程和接收线程结束
		recv_thread.join();//等待发送线程和接收线程结束
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