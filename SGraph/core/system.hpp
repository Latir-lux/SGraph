#ifndef SYSTEM_HPP
#define SYSTEM_HPP
#include <stdio.h>
#include <stdlib.h>
#include "core/engine.hpp"
#include "core/graph.hpp"
#include "core/filesystem.hpp"

template <typename VertexId, typename Weight>
struct UnionMessage
{
	std::pair<VertexId, Weight> msg[hubs];
};

template <typename Weight, typename EdgeData, typename Init, typename Forward, typename Backward, typename Merge>
class System
{
public:
	using Adjlist = std::vector<AdjEdge<EdgeData>>;
	Init init;
	Forward forward;
	Backward backward;
	Merge merge;

	uint32_t query_snapshot;

	bool symmetric;

	Graph<EdgeData> graph;
	Engine<EdgeData> query_engine;
	Engine<EdgeData> index_engine;
	MPI_Comm &index_comm;
	MPI_Comm &query_comm;

	std::vector<EdgeUnit<EdgeData>> add_buffer[2];
	std::vector<EdgeUnit<EdgeData>> del_buffer[2];

	// bad locality
	// VertexId (*incoming_indexes_parent[2])[hubs]; // VertexTree [2][vertices][hubs]
	// VertexId (*outgoing_indexes_parent[2])[hubs]; // VertexTree [2][vertices][hubs]
	// Weight (*incoming_indexes_value[2])[hubs];
	// Weight (*outgoing_indexes_value[2])[hubs];

	VertexId **incoming_indexes_parent[2]; // VertexTree [2][vertices][hubs]
	VertexId **outgoing_indexes_parent[2]; // VertexTree [2][vertices][hubs]
	Weight **incoming_indexes_value[2];
	Weight **outgoing_indexes_value[2];//在ss_query函数中计算

	std::mutex *index_mutex;

	VertexId *hub; // 记录所有中心节点，数量是hubs
	int *is_hub;   // 记录一个节点是否是中心节点，如果不是则为-1，如果是则为中心节点的编号（0~hubs-1）
	VertexSubset *active_all;
	VertexSubset *active_in;
	VertexSubset *active_out;
	VertexSubset *active_tree;
	VertexSubset *active_ins[hubs];
	VertexSubset *active_outs[hubs];
	VertexSubset *active_trees[hubs];

	System(MPI_Instance &mpi, VertexId _vertices, bool _symmetric = false, int _query_threads = 16, int _index_threads = 16)
		: query_comm(mpi.query_comm), index_comm(mpi.index_comm), graph(_vertices, _symmetric, query_snapshot), query_engine(&graph, _query_threads, QueryMessage, mpi.query_comm), index_engine(&graph, _index_threads, IndexMessage, mpi.index_comm)
	{
		omp_set_dynamic(0);//禁止动态调整线程数
		query_snapshot = 0;
		symmetric = _symmetric;
		hub = new VertexId[hubs];
		active_all = graph.alloc_vertex_subset();
		active_all->fill();
		active_in = graph.alloc_vertex_subset();
		active_out = graph.alloc_vertex_subset();
		active_tree = graph.alloc_vertex_subset();
		for (int i = 0; i < hubs; i++)
		{
			active_ins[i] = graph.alloc_vertex_subset();
			active_outs[i] = graph.alloc_vertex_subset();
			active_trees[i] = graph.alloc_vertex_subset();
		}
	}

	// used to evaluate larger graphs with limited memory
	// init load
	void load_file(std::string file)
	{
		long total_bytes = file_size(file.c_str(), graph.partition_id);
		size_t edge_unit_size = sizeof(EdgeUnit<EdgeData>);//文件中的边数total_edges是文件的字节大小除以每个边单位的大小。
		long total_edges = total_bytes / edge_unit_size;
		long vector_max_length = 1 << 26;// 2^26,long型总共32位,表示每次从文件中一次性读取的最大边的数量。
		long read_edges = 0;
		while (read_edges < total_edges)
		{//分块读取边并获取它们的度数，一个块的大小是2^26
			long current_read_edges = vector_max_length;
			if (current_read_edges > total_edges - read_edges)
			{
				current_read_edges = total_edges - read_edges;
			}
			EdgeUnit<EdgeData> *edges = new EdgeUnit<EdgeData>[current_read_edges];
			double t1 = -get_time();
			get_edge_vector(file, edges, read_edges, read_edges + current_read_edges, graph.partition_id);
			t1 += get_time();

			double t2 = -get_time();
			graph.init_get_degree(edges, current_read_edges);//读入的信息会存放到graph类中
			t2 += get_time();

			delete edges;
			read_edges += current_read_edges;
			if (graph.partition_id == 0)
			{
				printf("get degree: %ld / %ld\n", read_edges, total_edges);
				printf("t1 = %f t2 = %f\n", t1, t2);
			}
		}

		for (int v_i = 0; v_i < graph.vertices; v_i++)
		{
			graph.outgoing_storage[v_i]->adjlist.reserve(graph.out_degree[v_i]);
			graph.incoming_storage[v_i]->adjlist.reserve(graph.in_degree[v_i]);
		}
		//这段代码是为每个顶点的邻接边列表预留内存。它依赖于已经计算好的度数信息，以确保为邻接边列表预留的内存足够。这是为什么我们首先进行度数的计算的原因。
		//代码中对于边数据读取了两次，第一次读取的目的主要是计算每个顶点的度数（即每个顶点连接的边的数量）。第一次读取的目的主要是计算每个顶点的度数（即每个顶点连接的边的数量）。
		//这样做的好处是，当我们在第二次读取中实际地加载边到图数据结构中时，内存已经被预留，所以添加操作更加高效。如果没有这种预留，每次添加新的邻接边时，都可能需要动态地重新分配内存，这可能会导致大量的内存复制操作，从而影响性能。

		read_edges = 0;
		while (read_edges < total_edges)
		{
			long current_read_edges = vector_max_length;
			if (current_read_edges > total_edges - read_edges)
			{
				current_read_edges = total_edges - read_edges;
			}
			EdgeUnit<EdgeData> *edges = new EdgeUnit<EdgeData>[current_read_edges];
			double t1 = -get_time();
			get_edge_vector(file, edges, read_edges, read_edges + current_read_edges, graph.partition_id);
			t1 += get_time();

			double t2 = -get_time();
			graph.init_load_edges_visible(edges, current_read_edges);
			t2 += get_time();

			delete edges;
			read_edges += current_read_edges;
			if (graph.partition_id == 0)
			{
				printf("load edges: %ld / %ld\n", read_edges, total_edges);
				printf("t1 = %f t2 = %f\n", t1, t2);
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, graph.out_degree, graph.vertices, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, graph.in_degree, graph.vertices, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}

	// used to evaluate larger graphs with limited memory
	// init load (also prepare edges for update: 100*(5000add+5000del) 100*(50000add+50000del) 100*(500000add+500000del) 300 add_array and 300 del_array in total)
	void load_file_update(std::string file, EdgeUnit<EdgeData> **&dynamic_edges_add, EdgeUnit<EdgeData> **&dynamic_edges_del)
	{
		long total_bytes = file_size(file.c_str(), graph.partition_id);
		size_t edge_unit_size = sizeof(EdgeUnit<EdgeData>);
		long total_edges = total_bytes / edge_unit_size * 0.7;//使用文件大小除以 EdgeUnit<EdgeData> 的大小得到边的数量。这里只加载了文件中的70%边数据。
		long vector_max_length = 1 << 26;
		long read_edges = 0;
		while (read_edges < total_edges)
		{//在 while 循环中，文件中的边被分批加载。每次加载的边数量为 vector_max_length，或者如果剩余的边少于这个数量，则加载所有剩余的边。
			long current_read_edges = vector_max_length;
			if (current_read_edges > total_edges - read_edges)
			{
				current_read_edges = total_edges - read_edges;
			}
			EdgeUnit<EdgeData> *edges = new EdgeUnit<EdgeData>[current_read_edges];
			double t1 = -get_time();
			get_edge_vector(file, edges, read_edges, read_edges + current_read_edges, graph.partition_id);//get_edge_vector函数从文件中读取边。
			t1 += get_time();

			double t2 = -get_time();
			graph.init_get_degree(edges, current_read_edges);//对于每一批加载的边，它们的度数都会被初始化
			t2 += get_time();

			delete edges;
			read_edges += current_read_edges;
			if (graph.partition_id == 0)
			{
				printf("get degree: %ld / %ld\n", read_edges, total_edges);
				printf("t1 = %f t2 = %f\n", t1, t2);
			}
		}

		long dynamic_edge_num = 100 * (10000 + 100000 + 1000000);
		EdgeUnit<EdgeData> *dynamic_edges = new EdgeUnit<EdgeData>[dynamic_edge_num];
		get_edge_vector(file, dynamic_edges, total_edges, total_edges + dynamic_edge_num, graph.partition_id);
		graph.init_get_degree(dynamic_edges, dynamic_edge_num);

		for (int v_i = 0; v_i < graph.vertices; v_i++)
		{
			graph.outgoing_storage[v_i]->adjlist.reserve(graph.out_degree[v_i]);
			graph.incoming_storage[v_i]->adjlist.reserve(graph.in_degree[v_i]);
		}

		dynamic_edges_add = new EdgeUnit<EdgeData> *[300];
		dynamic_edges_del = new EdgeUnit<EdgeData> *[300];
		int cnt = 0;
		for (int i = 0; i < 300; i++)
		{
			int array_size;
			if (i / 100 == 0)
				array_size = 5000;
			if (i / 100 == 1)
				array_size = 50000;
			if (i / 100 == 2)
				array_size = 500000;
			dynamic_edges_add[i] = new EdgeUnit<EdgeData>[array_size];
			dynamic_edges_del[i] = new EdgeUnit<EdgeData>[array_size];
			for (int j = 0; j < array_size; j++)
			{
				dynamic_edges_add[i][j] = dynamic_edges[cnt];
				cnt++;
			}
			for (int j = 0; j < array_size; j++)
			{
				dynamic_edges_del[i][j] = dynamic_edges[cnt];
				cnt++;
			}
			graph.init_load_edges_invisible(dynamic_edges_add[i], array_size);
			graph.init_load_edges_visible(dynamic_edges_del[i], array_size);
		}

		read_edges = 0;
		while (read_edges < total_edges)
		{
			long current_read_edges = vector_max_length;
			if (current_read_edges > total_edges - read_edges)
			{
				current_read_edges = total_edges - read_edges;
			}
			EdgeUnit<EdgeData> *edges = new EdgeUnit<EdgeData>[current_read_edges];
			double t1 = -get_time();
			get_edge_vector(file, edges, read_edges, read_edges + current_read_edges, graph.partition_id);
			t1 += get_time();

			double t2 = -get_time();
			graph.init_load_edges_visible(edges, current_read_edges);
			t2 += get_time();

			delete edges;
			read_edges += current_read_edges;
			if (graph.partition_id == 0)
			{
				printf("load edges: %ld / %ld\n", read_edges, total_edges);
				printf("t1 = %f t2 = %f\n", t1, t2);
			}
		}

		MPI_Allreduce(MPI_IN_PLACE, graph.out_degree, graph.vertices, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, graph.in_degree, graph.vertices, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//使用MPI_Allreduce来同步所有MPI进程中的出度和入度信息，确保所有进程都有完整的度数信息。
	}

	// used to evaluate larger graphs with limited memory
	// update
	void sim_update(EdgeUnit<EdgeData> *dynamic_edges_add, EdgeUnit<EdgeData> *dynamic_edges_del, int array_size)
	{
		add_buffer[0].clear();
		add_buffer[1].clear();
		del_buffer[0].clear();
		del_buffer[1].clear();
		for (int i = 0; i < array_size; i++)
		{
			add_buffer[0].emplace_back(dynamic_edges_add[i]);
			add_buffer[1].emplace_back(dynamic_edges_add[i]);
			del_buffer[0].emplace_back(dynamic_edges_del[i]);
			del_buffer[1].emplace_back(dynamic_edges_del[i]);
		}
		graph.update_edges_invisible_to_visible(dynamic_edges_add, array_size);
		graph.update_edges_visible_to_invisible(dynamic_edges_del, array_size);
	}

	void add_edges(std::vector<EdgeUnit<EdgeData>> &edges)
	{
		for (auto it : edges)
		{
			add_buffer[0].emplace_back(it);
			add_buffer[1].emplace_back(it);
		}
		graph.add_edges(edges);
	}

	void del_edges(std::vector<EdgeUnit<EdgeData>> &edges)
	{
		for (auto it : edges)
		{
			del_buffer[0].emplace_back(it);
			del_buffer[1].emplace_back(it);
		}
		graph.del_edges(edges);
	}

	void step()
	{//system中的step函数会调用graph中的step函数，graph中的step函数会调用storage中的step函数，函数的作用是将buffer中的值更新到Adjlist中。
		graph.step();
		query_snapshot++;
	}

	inline bool update_index(VertexId &old_parent, Weight &old_value, VertexId new_parent, Weight new_value, std::mutex &mutex)
	{
		std::lock_guard<std::mutex> lock(mutex);
		if (merge(old_value, new_value) != old_value)
		{
			old_parent = new_parent;
			old_value = new_value;
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool monotonic_write(Weight *ptr, Weight val)
	{
		volatile Weight curr_val;
		bool done = false;
		do
		{
			curr_val = *ptr;
		} while ((merge(curr_val, val) != curr_val) && !(done = cas(ptr, curr_val, val)));
		return done;
	}

	void ss_query(int h_i, bool reverse = false)
	{//点对点查询，会获得所有顶点到中心点h_i的最短距离,以及最短路径上的父节点
		int i = (query_snapshot & 1) ^ 1;// `i`确定查询快照，可能为0或1
		VertexId **index_parent = reverse ? incoming_indexes_parent[i] : outgoing_indexes_parent[i];//index_parent存储了每个顶点到达每一个hub节点的路径上的父节点，根据是否逆转选择索引
		Weight **index_value = reverse ? incoming_indexes_value[i] : outgoing_indexes_value[i];//index_value是任意顶点到中心顶点的距离
		VertexId root = hub[h_i];// 设定root为指定的hub
		active_in->clear();
		active_in->set_bit(root);
		index_engine.template process_vertices<VertexId>(// 为每个顶点初始化其在索引中的值
			[&](VertexId vtx)
			{
				Weight init_val = reverse ? init(vtx, root) : init(root, vtx);
				index_parent[vtx][h_i] = graph.vertices;
				index_value[vtx][h_i] = init_val;
				return 1;
			},
			active_all);
		VertexId active_vertices = 1;
		for (int i_i = 0; active_vertices > 0; i_i++)
		{// 为每个顶点初始化其在索引中的值
			double start = get_time();
			active_out->clear();

			// test no buffer vs message buffer (when single node)
			// active_vertices = 0;
			// #pragma omp parallel for reduction(+:active_vertices)
			// for (VertexId v_i=0;v_i<graph.vertices;v_i++)
			// {
			// 	if (!active_in->get_bit(v_i)) continue;
			// 	VertexId local_active_vertices = 0;
			// 	uint32_t i = (query_snapshot & 1 ) ^ 1;
			// 	for (auto iter : index_engine.outgoing_adjlist[v_i]) {
			// 		if (!iter.get_valid(i)) continue;
			// 		VertexId dst = iter.nbr;
			// 		Weight relax_val = forward(index_value[v_i][h_i], iter.data);
			// 		if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i]) {
			// 			if (update_index(index_parent[dst][h_i], index_value[dst][h_i], v_i, relax_val, index_mutex[dst])) {
			// 				active_out->set_bit(dst);
			// 				local_active_vertices +=1;
			// 			}
			// 		}
			// 	}
			// 	active_vertices += local_active_vertices;
			// }

			//process_edge函数的定义：R process_edges(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, std::function<void(VertexId, Adjlist&)> dense_signal, std::function<R(VertexId, M)> dense_slot, Bitmap * active, Bitmap * dense_selective = nullptr, Bitmap ** dense_selective_hubs = nullptr) 
			active_vertices = index_engine.template process_edges<VertexId, std::pair<VertexId, Weight>>(
				[&](VertexId src)
				{//对于每个源顶点src，此回调发射/发送一个消息，该消息包含顶点自身及其与当前hub（由h_i表示）相关联的索引值。
					index_engine.emit(src, std::make_pair(src, index_value[src][h_i]));
				},
				[&](VertexId src, std::pair<VertexId, Weight> msg, Adjlist &outgoing_adj)
				{//这个回调对于每个源顶点src以及与它关联的消息msg和出边列表outgoing_adj，都会执行一些操作。它首先会计算一个值，称为relax_val，并尝试合并这个值到目的顶点dst的索引值。如果成功合并（即更新了索引值），它会标记dst为激活，并递增激活的计数。最后，这个回调返回激活的顶点数。
					VertexId activated = 0;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : outgoing_adj)
					{
						if (!iter.get_valid(i))
							continue;
						VertexId dst = iter.nbr;
						Weight relax_val = forward(msg.second, iter.data);
						if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i])
						{
							if (update_index(index_parent[dst][h_i], index_value[dst][h_i], msg.first, relax_val, index_mutex[dst]))
							{
								active_out->set_bit(dst);
								activated += 1;
							}
						}
					}
					return activated;
				},
				[&](VertexId dst, Adjlist &incoming_adj)
				{//这个回调对于每个源顶点src以及与它关联的消息msg和出边列表outgoing_adj，都会执行一些操作。它首先会计算一个值，称为relax_val，并尝试合并这个值到目的顶点dst的索引值。如果成功合并（即更新了索引值），它会标记dst为激活，并递增激活的计数。最后，这个回调返回激活的顶点数。
					Weight init_val = reverse ? init(dst, root) : init(root, dst);
					Weight msg = init_val;
					VertexId msg_src = graph.vertices;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : incoming_adj)
					{
						if (!iter.get_valid(i))
							continue;
						VertexId src = iter.nbr;
						// if (active_in->get_bit(src)) {
						Weight relax_val = forward(index_value[src][h_i], iter.data);
						if (merge(relax_val, msg) != msg)
						{
							msg = relax_val;
							msg_src = src;
						}
						// }
					}
					if (merge(msg, init_val) != init_val)
						index_engine.emit(dst, std::make_pair(msg_src, msg));
				},
				[&](VertexId dst, std::pair<VertexId, Weight> msg)
				{
					if (merge(msg.second, index_value[dst][h_i]) != index_value[dst][h_i])
					{
						index_parent[dst][h_i] = msg.first;
						index_value[dst][h_i] = msg.second;
						active_out->set_bit(dst);
						return 1;
					}
					return 0;
				},
				active_in);

			std::swap(active_in, active_out);
			double end = get_time();
		}
	}

	void build_index()
	{
		if (graph.partition_id == 0)//只在id为0的分区build_index
		{
			printf("build_index\n");
		}

		// bad locality
		// if (symmetric) {
		// 	for (int i=0;i<2;i++) {
		// 		incoming_indexes_parent[i] = outgoing_indexes_parent[i] = graph.template alloc_vertex_array_adhoc_perhub<VertexId>();
		// 		incoming_indexes_value[i] = outgoing_indexes_value[i] = graph.template alloc_vertex_array_adhoc_perhub<Weight>();
		// 	}
		// }
		// else {
		// 	for (int i=0;i<2;i++) {
		// 		incoming_indexes_parent[i] = graph.template alloc_vertex_array_adhoc_perhub<VertexId>();
		// 		incoming_indexes_value[i] = graph.template alloc_vertex_array_adhoc_perhub<Weight>();
		// 		outgoing_indexes_parent[i] = graph.template alloc_vertex_array_adhoc_perhub<VertexId>();
		// 		outgoing_indexes_value[i] = graph.template alloc_vertex_array_adhoc_perhub<Weight>();
		// 	}
		// }

		if (symmetric)
		{
			for (int i = 0; i < 2; i++)//两个版本
			{
				incoming_indexes_parent[i] = outgoing_indexes_parent[i] = new VertexId *[graph.vertices];
				incoming_indexes_value[i] = outgoing_indexes_value[i] = new Weight *[graph.vertices];
				for (VertexId v_i = graph.partition_offset[graph.partition_id]; v_i < graph.partition_offset[graph.partition_id + 1]; v_i++)
				{
					incoming_indexes_parent[i][v_i] = new VertexId[hubs];
				}
				for (VertexId v_i = graph.partition_offset[graph.partition_id]; v_i < graph.partition_offset[graph.partition_id + 1]; v_i++)
				{
					incoming_indexes_value[i][v_i] = new Weight[hubs];
				}
			}
		}
		else
		{
			for (int i = 0; i < 2; i++)//两个版本
			{
				incoming_indexes_parent[i] = new VertexId *[graph.vertices];
				incoming_indexes_value[i] = new Weight *[graph.vertices];
				outgoing_indexes_parent[i] = new VertexId *[graph.vertices];
				outgoing_indexes_value[i] = new Weight *[graph.vertices];
				for (VertexId v_i = graph.partition_offset[graph.partition_id]; v_i < graph.partition_offset[graph.partition_id + 1]; v_i++)
				{
					incoming_indexes_parent[i][v_i] = new VertexId[hubs];
					outgoing_indexes_parent[i][v_i] = new VertexId[hubs];
				}
				for (VertexId v_i = graph.partition_offset[graph.partition_id]; v_i < graph.partition_offset[graph.partition_id + 1]; v_i++)
				{
					incoming_indexes_value[i][v_i] = new Weight[hubs];
					outgoing_indexes_value[i][v_i] = new Weight[hubs];
				}
			}
		}

		index_mutex = new std::mutex[graph.vertices];
		is_hub = graph.template alloc_vertex_array<int>();
		// select hubs
		VertexId *sum_degree = graph.template alloc_vertex_array<VertexId>();
#pragma omp parallel for
		for (VertexId v_i = 0; v_i < graph.vertices; v_i++)
		{
			sum_degree[v_i] = graph.out_degree[v_i] + graph.in_degree[v_i];
		};
#pragma omp parallel for
		for (VertexId v_i = 0; v_i < graph.vertices; v_i++)
		{
			is_hub[v_i] = -1;
		}
		for (int h_i = 0; h_i < hubs; h_i++)
		{//选择度数最大的hub个顶点，作为中心顶点
			VertexId max_v_i = graph.vertices;
			for (VertexId v_i = 0; v_i < graph.vertices; v_i++)
			{
				if (is_hub[v_i] != -1)
					continue;
				if (max_v_i == graph.vertices || sum_degree[v_i] > sum_degree[max_v_i])
				{
					max_v_i = v_i;
				}
			}
			is_hub[max_v_i] = h_i;
			hub[h_i] = max_v_i;
		}

		// cal hub's sssp
		double tot_time = 0;
		for (int h_i = 0; h_i < hubs; h_i++)
		{//统计所有顶点到hub的最短路径
			double t = -get_time();
			ss_query(h_i);////点对点查询，会获得所有顶点到中心点h_i的最短距离,以及最短路径上的父节点
			t += get_time();
			tot_time += t;
			if (graph.partition_id == 0)
			{
				printf("t = %f\n", t);
			}

			if (!symmetric)
			{
				double t = -get_time();
				index_engine.transpose();
				ss_query(h_i, true);//ss_query(h_i)执行的是false方向的查询，如果是对于非对称图（有向图），还需要在true方向再查一次。
				index_engine.transpose();
				t += get_time();
				tot_time += t;
				if (graph.partition_id == 0)
				{
					printf("t = %f\n", t);
				}
			}
		}

		int i = query_snapshot & 1;
#pragma omp parallel for
		for (VertexId v_i = graph.partition_offset[graph.partition_id]; v_i < graph.partition_offset[graph.partition_id + 1]; v_i++)
		{
			for (int h_i = 0; h_i < hubs; h_i++)
			{
				incoming_indexes_parent[i][v_i][h_i] = incoming_indexes_parent[i ^ 1][v_i][h_i];
				incoming_indexes_value[i][v_i][h_i] = incoming_indexes_value[i ^ 1][v_i][h_i];
				if (!symmetric)
				{
					outgoing_indexes_parent[i][v_i][h_i] = outgoing_indexes_parent[i ^ 1][v_i][h_i];
					outgoing_indexes_value[i][v_i][h_i] = outgoing_indexes_value[i ^ 1][v_i][h_i];
				}
			}
		}

		// printf("build_index ok\n");
	}

	/*
	void modify_index_add(int h_i, std::vector<EdgeUnit<EdgeData>>& edges, bool reverse = false){
		int i = (query_snapshot & 1) ^ 1;
		VertexId ** index_parent = reverse ? incoming_indexes_parent[i] : outgoing_indexes_parent[i];
		Weight ** index_value = reverse ? incoming_indexes_value[i] : outgoing_indexes_value[i];
		VertexId root = hub[h_i];

		VertexId active_vertices = 0;
		Weight * src_val = new Weight [edges.size()];
		#pragma omp parallel for
		for (int i=0;i<edges.size();i++) {
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (reverse) std::swap(src, dst);
			int src_part = graph.get_partition_id(src);
			if (graph.partition_id == src_part) {
				src_val[i] = index_value[src][h_i];
			}
			else {
				src_val[i] = 0;
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, src_val, edges.size(), get_mpi_data_type<Weight>(), MPI_SUM, index_comm);
		#pragma omp parallel for
		for (int i=0;i<edges.size();i++) {
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (reverse) std::swap(src, dst);
			int dst_part = graph.get_partition_id(dst);
			if (graph.partition_id == dst_part) {
				Weight edge_data = edges[i].edge_data;
				Weight relax_val = forward(edge_data, src_val[i]);
				if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i]) {
					if (update_index(index_parent[dst][h_i], index_value[dst][h_i], src, relax_val, index_mutex[dst])) {
						active_in->set_bit(dst);
						active_vertices++;
					}
				}
			}
		}
		delete [] src_val;
		MPI_Allreduce(MPI_IN_PLACE, &active_vertices, 1, MPI_INT, MPI_SUM, index_comm);


		for (int i_i=0;active_vertices>0;i_i++) {
			active_out->clear();
			active_vertices = index_engine.template process_edges<VertexId,std::pair<VertexId, Weight>>(
				[&](VertexId src){
					index_engine.emit(src, std::make_pair(src, index_value[src][h_i]));
				},
				[&](VertexId src, std::pair<VertexId, Weight> msg, Adjlist& outgoing_adj){
					VertexId activated = 0;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : outgoing_adj){
						if (!iter.get_valid(i)) continue;
						VertexId dst = iter.nbr;
						Weight relax_val = forward(msg.second, iter.data);
						if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i]) {
							if (update_index(index_parent[dst][h_i], index_value[dst][h_i], src, relax_val, index_mutex[dst])) {
								active_out->set_bit(dst);
								activated +=1;
							}
						}
					}
					return activated;
				},
				[&](VertexId dst, Adjlist& incoming_adj) {
					Weight init_val = reverse ? init(dst, root) : init(root, dst);
					Weight msg = init_val;
					VertexId msg_src = graph.vertices;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : incoming_adj) {
						if (!iter.get_valid(i)) continue;
						VertexId src = iter.nbr;
						// if (active_in->get_bit(src)) {
							Weight relax_val = forward(index_value[src][h_i], iter.data);
							if (merge(relax_val, msg) != msg) {
								msg = relax_val;
								msg_src = src;
							}
						// }
					}
					if (merge(msg, init_val) != init_val) index_engine.emit(dst, std::make_pair(msg_src, msg));
				},
				[&](VertexId dst, std::pair<VertexId, Weight>msg) {
					if (merge(msg.second, index_value[dst][h_i]) != index_value[dst][h_i]) {
						index_parent[dst][h_i] = msg.first;
						index_value[dst][h_i] = msg.second;
						active_out->set_bit(dst);
						return 1;
					}
					return 0;
				},
				active_in
			);
			std::swap(active_in, active_out);
		}
	}

	void modify_index_add(){
		uint32_t i = (query_snapshot & 1) ^ 1;
		if (add_buffer[i].size() == 0) return;
		for (int h_i=0;h_i<hubs;h_i++){
			modify_index_add(h_i, add_buffer[i]);
			if (!symmetric) {
				index_engine.transpose();
				modify_index_add(h_i, add_buffer[i], true);
				index_engine.transpose();
			}
		}
		add_buffer[i].clear();
	}

	void modify_index_del(int h_i, std::vector<EdgeUnit<EdgeData>>& edges, bool reverse = false){
		int i = (query_snapshot & 1) ^ 1;
		VertexId ** index_parent = reverse ? incoming_indexes_parent[i] : outgoing_indexes_parent[i];
		Weight ** index_value = reverse ? incoming_indexes_value[i] : outgoing_indexes_value[i];
		VertexId root = hub[h_i];

		// figure out invalid vertices
		active_in->clear();
		active_tree->clear();
		VertexId active_vertices = 0;
		// double time = -get_time();
		#pragma omp parallel for
		for (int i=0;i<edges.size();i++){
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (reverse) std::swap(src, dst);
			if (!(graph.partition_offset[graph.partition_id] <= dst && dst < graph.partition_offset[graph.partition_id + 1])) continue;
			if (index_parent[dst][h_i] == src){
				active_in->set_bit(dst);
				active_tree->set_bit(dst);
				active_vertices++;
			}
		}
		// time += get_time();
		// printf("prepare %f\n",time);
		// time = -get_time();
		MPI_Allreduce(MPI_IN_PLACE, &active_vertices, 1, MPI_INT, MPI_SUM, index_comm);

		if (graph.partition_id==0) {
			//printf("active_vertices = %d\n",active_vertices);
		}

		for (int i_i=0;active_vertices>0;i_i++){
			active_out->clear();
			active_vertices = index_engine.template process_edges_sparse<VertexId, VertexId>(
				[&](VertexId src){
					index_engine.emit(src, src);
				},
				[&](VertexId src, VertexId msg, Adjlist& outgoing_adj){
					VertexId activated = 0;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : outgoing_adj) {
						if (!iter.get_valid(i)) continue;
						VertexId dst = iter.nbr;
						if (index_parent[dst][h_i] == src){
							active_out->set_bit(dst);
							active_tree->set_bit(dst);
							activated++;
						}
					}
					return activated;
				},
				active_in
			);
			std::swap(active_in, active_out);
		}

		// initialize invalid vertices
		active_vertices = index_engine.template process_vertices<VertexId>(
			[&](VertexId vtx){
				Weight init_val = reverse ? init(vtx, root) : init(root, vtx);
				index_parent[vtx][h_i] = graph.vertices;
				index_value[vtx][h_i] = init_val;
				return 1;
			},
			active_tree
		);

		if (graph.partition_id==0) {
			//printf("active_vertices = %d\n",active_vertices);
		}

		// time += get_time();
		// printf("init %f\n",time);
		// time = -get_time();
		//if (active_vertices) printf("              active  %d\n", active_vertices);

		// recompute invalid vertices
		if (active_vertices){
			active_out->clear();
			active_vertices = index_engine.template process_edges_dense<VertexId,std::pair<VertexId, Weight>>(
				[&](VertexId dst, Adjlist& incoming_adj) {
					if (!active_tree->get_bit(dst)) return;
					Weight init_val = reverse ? init(dst, root) : init(root, dst);
					Weight msg = init_val;
					VertexId msg_src = graph.vertices;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : incoming_adj) {
						if (!iter.get_valid(i)) continue;
						VertexId src = iter.nbr;
						Weight relax_val = forward(index_value[src][h_i], iter.data);
						if (merge(relax_val, msg) != msg) {
							msg = relax_val;
							msg_src = src;
						}
					}
					if (msg_src != graph.vertices) index_engine.emit(dst, std::make_pair(msg_src, msg));
				},
				[&](VertexId dst, std::pair<VertexId, Weight> msg) {
					if (merge(msg.second, index_value[dst][h_i]) != index_value[dst][h_i]) {
						index_parent[dst][h_i] = msg.first;
						index_value[dst][h_i] = msg.second;
						active_out->set_bit(dst);
						return 1;
					}
					return 0;
				},
				active_all,
				active_tree
			);
			std::swap(active_in, active_out);
		}

		// time += get_time();
		// printf("pull %f\n",time);
		// time = -get_time();

		for (int i_i=0;active_vertices>0;i_i++) {
			active_out->clear();
			active_vertices = index_engine.template process_edges<VertexId,std::pair<VertexId, Weight>>(
				[&](VertexId src){
					index_engine.emit(src, std::make_pair(src, index_value[src][h_i]));
				},
				[&](VertexId src, std::pair<VertexId, Weight> msg, Adjlist& outgoing_adj){
					VertexId activated = 0;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : outgoing_adj) {
						if (!iter.get_valid(i)) continue;
						VertexId dst = iter.nbr;
						Weight relax_val = forward(msg.second, iter.data);
						if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i]) {
							if (update_index(index_parent[dst][h_i], index_value[dst][h_i], src, relax_val, index_mutex[dst])) {
								active_out->set_bit(dst);
								activated +=1;
							}
						}
					}
					return activated;
				},
				[&](VertexId dst, Adjlist& incoming_adj) {
					Weight init_val = reverse ? init(dst, root) : init(root, dst);
					Weight msg = init_val;
					VertexId msg_src = graph.vertices;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (auto iter : incoming_adj) {
						if (!iter.get_valid(i)) continue;
						VertexId src = iter.nbr;
						// if (active_in->get_bit(src)) {
							Weight relax_val = forward(index_value[src][h_i], iter.data);
							if (merge(relax_val, msg) != msg) {
								msg = relax_val;
								msg_src = src;
							}
						// }
					}
					if (merge(msg, init_val) != init_val) index_engine.emit(dst, std::make_pair(msg_src, msg));
				},
				[&](VertexId dst, std::pair<VertexId, Weight>msg) {
					if (merge(msg.second, index_value[dst][h_i]) != index_value[dst][h_i]) {
						index_parent[dst][h_i] = msg.first;
						index_value[dst][h_i] = msg.second;
						active_out->set_bit(dst);
						return 1;
					}
					return 0;
				},
				active_in
			);
			std::swap(active_in, active_out);
		}
		// time += get_time();
		// printf("interate %f\n",time);
	}

	void modify_index_del(){
		uint32_t i = (query_snapshot & 1) ^ 1;
		if (del_buffer[i].size() == 0) return;
		for (int h_i=0;h_i<hubs;h_i++){
			modify_index_del(h_i, del_buffer[i]);
			if (!symmetric) {
				index_engine.transpose();
				modify_index_del(h_i, del_buffer[i], true);
				index_engine.transpose();
			}
		}
		del_buffer[i].clear();
	}
	*/

	void modify_index_hybrid(std::vector<EdgeUnit<EdgeData>> &add_edges, std::vector<EdgeUnit<EdgeData>> &del_edges, bool reverse = false)
	{
		int i = (query_snapshot & 1) ^ 1;
		VertexId **index_parent = reverse ? incoming_indexes_parent[i] : outgoing_indexes_parent[i];
		Weight **index_value = reverse ? incoming_indexes_value[i] : outgoing_indexes_value[i];
		VertexId *root = hub;

		// figure out invalid vertices
		active_in->clear();
		active_tree->clear();
		for (int i = 0; i < hubs; i++)
		{
			active_ins[i]->clear();
			active_trees[i]->clear();
		}
		VertexId active_vertices = 0;
// double time = -get_time();
#pragma omp parallel for
		for (int i = 0; i < del_edges.size(); i++)
		{
			VertexId src = del_edges[i].src;
			VertexId dst = del_edges[i].dst;
			if (reverse)
				std::swap(src, dst);
			if (!(graph.partition_offset[graph.partition_id] <= dst && dst < graph.partition_offset[graph.partition_id + 1]))
				continue;
			for (int h_i = 0; h_i < hubs; h_i++)
			{
				if (index_parent[dst][h_i] == src)
				{
					active_ins[h_i]->set_bit(dst);
					active_trees[h_i]->set_bit(dst);
					active_tree->set_bit(dst);
					active_vertices++;
				}
			}
		}

		// time += get_time();
		// printf("prepare %f\n",time);
		// time = -get_time();

		for (int h_i = 0; h_i < hubs; h_i++)
		{
			for (int i_i = 0; active_vertices > 0; i_i++)
			{
				active_out->clear();
				active_vertices = index_engine.template process_edges_sparse<VertexId, VertexId>(
					[&](VertexId src)
					{
						index_engine.emit(src, src);
					},
					[&](VertexId src, VertexId msg, Adjlist &outgoing_adj)
					{
						VertexId activated = 0;
						uint32_t i = (query_snapshot & 1) ^ 1;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							if (!active_tree->get_bit(dst) && index_parent[dst][h_i] == src)
							{
								active_trees[h_i]->set_bit(dst);
								active_out->set_bit(dst);
								active_tree->set_bit(dst);
								activated++;
							}
						}
						return activated;
					},
					active_ins[h_i]);
				std::swap(active_ins[h_i], active_out);
			}
		}

		// initialize invalid vertices
		active_vertices += index_engine.template process_vertices<VertexId>(
			[&](VertexId vtx)
			{
				for (int h_i = 0; h_i < hubs; h_i++)
				{
					if (active_trees[h_i]->get_bit(vtx))
					{
						Weight init_val = reverse ? init(vtx, root[h_i]) : init(root[h_i], vtx);
						index_parent[vtx][h_i] = graph.vertices;
						index_value[vtx][h_i] = init_val;
					}
				}
				return 1;
			},
			active_tree);

		if (graph.partition_id == 0)
		{
			// printf("active_vertices = %d\n",active_vertices);
		}

		// time += get_time();
		// printf("init %f\n",time);
		// time = -get_time();
		// if (active_vertices) printf("              active  %d\n", active_vertices);

		MPI_Barrier(MPI_COMM_WORLD);
		double t = -get_time();

		// recompute invalid vertices
		if (active_vertices)
		{
			active_out->clear();
			for (int i = 0; i < hubs; i++)
			{
				active_outs[i]->clear();
			}
			active_vertices = index_engine.template process_edges_dense<VertexId, UnionMessage<VertexId, Weight>>(
				[&](VertexId dst, Adjlist &incoming_adj)
				{
					if (!active_tree->get_bit(dst))
						return;
					UnionMessage<VertexId, Weight> msgs;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						if (!active_trees[h_i]->get_bit(dst))
						{
							msgs.msg[h_i] = std::make_pair(graph.vertices, -1);
							continue;
						}
						Weight init_val = reverse ? init(dst, root[h_i]) : init(root[h_i], dst);
						Weight msg = init_val;
						VertexId msg_src = graph.vertices;
						uint32_t i = (query_snapshot & 1) ^ 1;
						for (auto iter : incoming_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId src = iter.nbr;
							// if (active_in->get_bit(src)) {
							Weight relax_val = forward(index_value[src][h_i], iter.data);
							if (merge(relax_val, msg) != msg)
							{
								msg = relax_val;
								msg_src = src;
							}
							// }
						}
						if (msg_src != graph.vertices)
						{
							msgs.msg[h_i] = std::make_pair(msg_src, msg);
						}
						else
						{
							msgs.msg[h_i] = std::make_pair(graph.vertices, -1);
						}
					}
					index_engine.emit(dst, msgs);
				},
				[&](VertexId dst, UnionMessage<VertexId, Weight> msgs)
				{
					VertexId activated = 0;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						if (msgs.msg[h_i].first == graph.vertices)
							continue;
						if (merge(msgs.msg[h_i].second, index_value[dst][h_i]) != index_value[dst][h_i])
						{
							index_parent[dst][h_i] = msgs.msg[h_i].first;
							index_value[dst][h_i] = msgs.msg[h_i].second;
							active_out->set_bit(dst);
							active_outs[h_i]->set_bit(dst);
							activated++;
						}
					}
					return activated;
				},
				active_all,
				active_tree,
				active_trees);
			std::swap(active_in, active_out);
			for (int i = 0; i < hubs; i++)
			{
				std::swap(active_ins[i], active_outs[i]);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		t += get_time();
		// if (graph.partition_id==0) {
		// 	printf("pull %f\n", t);
		// }

		// time += get_time();
		// printf("pull %f\n",time);
		// time = -get_time();

		if (add_edges.size() != 0)
		{
			VertexId active_vertices1 = 0;
			Weight *src_val[hubs];
			for (int i = 0; i < hubs; i++)
			{
				src_val[i] = new Weight[add_edges.size()];
			}
#pragma omp parallel for
			for (int i = 0; i < add_edges.size(); i++)
			{
				VertexId src = add_edges[i].src;
				VertexId dst = add_edges[i].dst;
				if (reverse)
					std::swap(src, dst);
				int src_part = graph.get_partition_id(src);
				if (graph.partition_id == src_part)
				{
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						src_val[h_i][i] = index_value[src][h_i];
					}
				}
				else
				{
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						src_val[h_i][i] = 0;
					}
				}
			}
			for (int i = 0; i < hubs; i++)
			{
				MPI_Allreduce(MPI_IN_PLACE, src_val[i], add_edges.size(), get_mpi_data_type<Weight>(), MPI_SUM, index_comm);
			}
#pragma omp parallel for
			for (int i = 0; i < add_edges.size(); i++)
			{
				VertexId src = add_edges[i].src;
				VertexId dst = add_edges[i].dst;
				if (reverse)
					std::swap(src, dst);
				int dst_part = graph.get_partition_id(dst);
				if (graph.partition_id == dst_part)
				{
					Weight edge_data = add_edges[i].edge_data;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						Weight relax_val = forward(edge_data, src_val[h_i][i]);
						if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i])
						{
							if (update_index(index_parent[dst][h_i], index_value[dst][h_i], src, relax_val, index_mutex[dst]))
							{
								active_in->set_bit(dst);
								active_ins[h_i]->set_bit(dst);
								active_vertices1++;
							}
						}
					}
				}
			}
			for (int i = 0; i < hubs; i++)
			{
				delete[] src_val[i];
			}
			MPI_Allreduce(MPI_IN_PLACE, &active_vertices1, 1, MPI_INT, MPI_SUM, index_comm);
			active_vertices += active_vertices1;
		}

		for (int i_i = 0; active_vertices > 0; i_i++)
		{
			active_out->clear();
			for (int i = 0; i < hubs; i++)
			{
				active_outs[i]->clear();
			}
			active_vertices = index_engine.template process_edges<VertexId, UnionMessage<VertexId, Weight>>(
				[&](VertexId src)
				{
					UnionMessage<VertexId, Weight> msgs;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						if (active_ins[h_i]->get_bit(src))
							msgs.msg[h_i] = std::make_pair(src, index_value[src][h_i]);
						else
							msgs.msg[h_i] = std::make_pair(graph.vertices, -1);
					}
					index_engine.emit(src, msgs);
				},
				[&](VertexId src, UnionMessage<VertexId, Weight> msgs, Adjlist &outgoing_adj)
				{
					VertexId activated = 0;
					uint32_t i = (query_snapshot & 1) ^ 1;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						if (msgs.msg[h_i].first == graph.vertices)
							continue;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							Weight relax_val = forward(msgs.msg[h_i].second, iter.data);
							if (merge(relax_val, index_value[dst][h_i]) != index_value[dst][h_i])
							{
								if (update_index(index_parent[dst][h_i], index_value[dst][h_i], src, relax_val, index_mutex[dst]))
								{
									active_out->set_bit(dst);
									active_outs[h_i]->set_bit(dst);
									activated += 1;
								}
							}
						}
					}
					return activated;
				},
				[&](VertexId dst, Adjlist &incoming_adj)
				{
					UnionMessage<VertexId, Weight> msgs;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						Weight init_val = reverse ? init(dst, root[h_i]) : init(root[h_i], dst);
						Weight msg = init_val;
						VertexId msg_src = graph.vertices;
						uint32_t i = (query_snapshot & 1) ^ 1;
						for (auto iter : incoming_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId src = iter.nbr;
							// if (active_in->get_bit(src)) {
							Weight relax_val = forward(index_value[src][h_i], iter.data);
							if (merge(relax_val, msg) != msg)
							{
								msg = relax_val;
								msg_src = src;
							}
							// }
						}
						if (merge(msg, init_val) != init_val)
							msgs.msg[h_i] = std::make_pair(msg_src, msg);
						else
							msgs.msg[h_i] = std::make_pair(graph.vertices, -1);
					}
					index_engine.emit(dst, msgs);
				},
				[&](VertexId dst, UnionMessage<VertexId, Weight> msgs)
				{
					VertexId activated = 0;
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						if (msgs.msg[h_i].first == graph.vertices)
							continue;
						if (merge(msgs.msg[h_i].second, index_value[dst][h_i]) != index_value[dst][h_i])
						{
							index_parent[dst][h_i] = msgs.msg[h_i].first;
							index_value[dst][h_i] = msgs.msg[h_i].second;
							active_out->set_bit(dst);
							active_outs[h_i]->set_bit(dst);
							activated++;
						}
					}
					return activated;
				},
				active_in);
			std::swap(active_in, active_out);
			for (int i = 0; i < hubs; i++)
			{
				std::swap(active_ins[i], active_outs[i]);
			}
		}

		// time += get_time();
		// printf("interate %f\n",time);
	}

	void modify_index()
	{
		uint32_t i = (query_snapshot & 1) ^ 1;
		modify_index_hybrid(add_buffer[i], del_buffer[i]);
		if (!symmetric)
		{
			index_engine.transpose();
			modify_index_hybrid(add_buffer[i], del_buffer[i], true);
			index_engine.transpose();
		}
		add_buffer[i].clear();
		del_buffer[i].clear();
	}

	Weight pnp(VertexId source, VertexId sink, uint64_t &active_num, uint64_t &iteration_num, uint64_t *active_on_iteration)
	{
		int i = query_snapshot & 1;
		Weight *value_source = graph.template alloc_vertex_array_adhoc<Weight>();
		Weight *value_sink = graph.template alloc_vertex_array_adhoc<Weight>();
		VertexSubset *active_source_in = graph.alloc_vertex_subset();
		VertexSubset *active_source_out = graph.alloc_vertex_subset();
		VertexSubset *active_sink_in = graph.alloc_vertex_subset();
		VertexSubset *active_sink_out = graph.alloc_vertex_subset();
		active_source_in->clear();
		active_source_in->set_bit(source);
		active_sink_in->clear();
		active_sink_in->set_bit(sink);

		query_engine.template process_vertices<VertexId>(
			[&](VertexId vtx)
			{
				value_source[vtx] = init(source, vtx);
				value_sink[vtx] = init(vtx, sink);
				return 0;
			},
			active_all);
		VertexId active_source_vertices = 1;
		VertexId active_sink_vertices = 1;
		Weight upper_bound = init(source, sink);
		int direction = 0; // 0:bidirection 1:forward 2:backward
		int intersect = 0;

		for (int i_i = 0; active_source_vertices > 0 || active_sink_vertices > 0; i_i++)
		{
			if (active_source_vertices > 0 && direction != 2)
			{ // forward
				iteration_num++;
				active_num += active_source_vertices;
				active_on_iteration[i_i] += active_source_vertices;
				active_source_out->clear();
				//R process_edges(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, std::function<void(VertexId, Adjlist&)> dense_signal, std::function<R(VertexId, M)> dense_slot, Bitmap * active, Bitmap * dense_selective = nullptr, Bitmap ** dense_selective_hubs = nullptr) {
				active_source_vertices = query_engine.template process_edges<VertexId, Weight>(
					[&](VertexId src)
					{
						query_engine.emit(src, value_source[src]);
					},
					[&](VertexId src, Weight msg, Adjlist &outgoing_adj)
					{
						VertexId activated = 0;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							Weight relax_val = forward(msg, iter.data);
							if (merge(relax_val, value_source[dst]) != value_source[dst])
							{
								if (monotonic_write(&value_source[dst], relax_val))
								{
									Weight tmp = forward(value_source[dst], value_sink[dst]);
									if (merge(tmp, init(source, sink)) != init(source, sink))
									{
										intersect = 1;
									}
									if (merge(tmp, upper_bound) != upper_bound)
									{
										monotonic_write(&upper_bound, tmp);
									}
									if (merge(value_source[dst], upper_bound) != upper_bound)
									{
										active_source_out->set_bit(dst);
										activated += 1;
									}
								}
							}
						}
						return activated;
					},
					[&](VertexId dst, Adjlist &incoming_adj)
					{
						Weight msg = init(source, sink);
						for (auto iter : incoming_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId src = iter.nbr;
							if (active_source_in->get_bit(src))
							{
								Weight relax_val = forward(value_source[src], iter.data);
								if (merge(relax_val, msg) != msg)
								{
									msg = relax_val;
								}
							}
						}
						if (merge(msg, upper_bound) != upper_bound)
							query_engine.emit(dst, msg); // trim by current shortest path
					},
					[&](VertexId dst, Weight msg)
					{
						if (merge(msg, value_source[dst]) != value_source[dst])
						{
							monotonic_write(&value_source[dst], msg);
							Weight tmp = forward(value_source[dst], value_sink[dst]);
							if (merge(tmp, init(source, sink)) != init(source, sink))
							{
								intersect = 1;
							}
							if (merge(tmp, upper_bound) != upper_bound)
							{
								monotonic_write(&upper_bound, tmp);
							}
							if (merge(value_source[dst], upper_bound) != upper_bound)
							{
								active_source_out->set_bit(dst);
								return 1;
							}
						}
						return 0;
					},
					active_source_in);
				std::swap(active_source_in, active_source_out);
				active_source_vertices = query_engine.template process_vertices<VertexId>(
					[&](VertexId vtx)
					{
						if (active_source_in->get_bit(vtx))
							return 1;
						return 0;
					},
					active_all);
#ifdef BFS_OPT
				MPI_Allreduce(MPI_IN_PLACE, &intersect, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
				if (intersect)
				{
					break;
				}
#endif
			}

			if (active_sink_vertices > 0 && direction != 1)
			{ // backward
				iteration_num++;
				active_num += active_sink_vertices;
				active_on_iteration[i_i] += active_sink_vertices;
				query_engine.transpose();
				active_sink_out->clear();
				active_sink_vertices = query_engine.template process_edges<VertexId, Weight>(
					[&](VertexId src)
					{
						query_engine.emit(src, value_sink[src]);
					},
					[&](VertexId src, Weight msg, Adjlist &outgoing_adj)
					{
						VertexId activated = 0;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							Weight relax_val = forward(msg, iter.data);
							if (merge(relax_val, value_sink[dst]) != value_sink[dst])
							{
								if (monotonic_write(&value_sink[dst], relax_val))
								{
									Weight tmp = forward(value_source[dst], value_sink[dst]);
									if (merge(tmp, init(source, sink)) != init(source, sink))
									{
										intersect = 1;
									}
									if (merge(tmp, upper_bound) != upper_bound)
									{
										monotonic_write(&upper_bound, tmp);
									}
									if (merge(value_sink[dst], upper_bound) != upper_bound)
									{
										active_sink_out->set_bit(dst);
										activated += 1;
									}
								}
							}
						}
						return activated;
					},
					[&](VertexId dst, Adjlist &incoming_adj)
					{
						Weight msg = init(source, sink);
						for (auto iter : incoming_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId src = iter.nbr;
							if (active_sink_in->get_bit(src))
							{
								Weight relax_val = forward(value_sink[src], iter.data);
								if (merge(relax_val, msg) != msg)
								{
									msg = relax_val;
								}
							}
						}
						if (merge(msg, upper_bound) != upper_bound)
							query_engine.emit(dst, msg); // trim by current shortest path
					},
					[&](VertexId dst, Weight msg)
					{
						if (merge(msg, value_sink[dst]) != value_sink[dst])
						{
							monotonic_write(&value_sink[dst], msg);
							Weight tmp = forward(value_source[dst], value_sink[dst]);
							if (merge(tmp, init(source, sink)) != init(source, sink))
							{
								intersect = 1;
							}
							if (merge(tmp, upper_bound) != upper_bound)
							{
								monotonic_write(&upper_bound, tmp);
							}
							if (merge(value_sink[dst], upper_bound) != upper_bound)
							{
								active_sink_out->set_bit(dst);
								return 1;
							}
						}
						return 0;
					},
					active_sink_in);
				query_engine.transpose();
				std::swap(active_sink_in, active_sink_out);
				active_sink_vertices = query_engine.template process_vertices<VertexId>(
					[&](VertexId vtx)
					{
						if (active_sink_in->get_bit(vtx))
							return 1;
						return 0;
					},
					active_all);
#ifdef BFS_OPT
				MPI_Allreduce(MPI_IN_PLACE, &intersect, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
				if (intersect)
				{
					break;
				}
#endif
			}

			// share upper bound
			if (graph.partition_id == 0)
			{
				for (int p_i = 1; p_i < graph.partitions; p_i++)
				{
					Weight upper_bound_remote;
					MPI_Status status;
					MPI_Recv(&upper_bound_remote, 1, get_mpi_data_type<Weight>(), p_i, QueryMessage, query_comm, &status);
					upper_bound = merge(upper_bound, upper_bound_remote);
				}
			}
			else
			{
				MPI_Send(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, QueryMessage, query_comm);
			}
			MPI_Bcast(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, query_comm);

			MPI_Allreduce(MPI_IN_PLACE, &intersect, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			if (direction == 0 && intersect)
			{
				if (active_source_vertices < active_sink_vertices)
				{
					direction = 1;
					active_sink_vertices = 0;
				}
				else
				{
					direction = 2;
					active_source_vertices = 0;
				}
			}
			if (direction == 0 && (active_source_vertices == 0 || active_sink_vertices == 0))
			{
				break;
			}
		}
		graph.dealloc_vertex_array_adhoc(value_source);
		graph.dealloc_vertex_array_adhoc(value_sink);
		delete active_source_in;
		delete active_source_out;
		delete active_sink_in;
		delete active_sink_out;
		return upper_bound;
	}

	Weight tripoline(VertexId source, VertexId sink, uint64_t &active_num, uint64_t &iteration_num, uint64_t *active_on_iteration)
	{
		int i = query_snapshot & 1;
		Weight **outgoing_index_value = outgoing_indexes_value[i];
		Weight **incoming_index_value = incoming_indexes_value[i];

		Weight *outgoing_index_value_source = new Weight[hubs];
		Weight *incoming_index_value_source = new Weight[hubs];
		Weight *outgoing_index_value_sink = new Weight[hubs];
		Weight *incoming_index_value_sink = new Weight[hubs];
		int source_part = graph.get_partition_id(source);
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			if (source_part == graph.partition_id)
			{
				outgoing_index_value_source[h_i] = outgoing_index_value[source][h_i];
				incoming_index_value_source[h_i] = incoming_index_value[source][h_i];
			}
		}
		MPI_Bcast(outgoing_index_value_source, hubs * sizeof(Weight), MPI_CHAR, source_part, query_comm);
		MPI_Bcast(incoming_index_value_source, hubs * sizeof(Weight), MPI_CHAR, source_part, query_comm);

		int sink_part = graph.get_partition_id(sink);
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			if (sink_part == graph.partition_id)
			{
				outgoing_index_value_sink[h_i] = outgoing_index_value[sink][h_i];
				incoming_index_value_sink[h_i] = incoming_index_value[sink][h_i];
			}
		}
		MPI_Bcast(outgoing_index_value_sink, hubs * sizeof(Weight), MPI_CHAR, sink_part, query_comm);
		MPI_Bcast(incoming_index_value_sink, hubs * sizeof(Weight), MPI_CHAR, sink_part, query_comm);

		Weight *value_source = graph.template alloc_vertex_array_adhoc<Weight>();
		VertexSubset *active_source_in = graph.alloc_vertex_subset();
		VertexSubset *active_source_out = graph.alloc_vertex_subset();
		active_source_in->clear();
		active_source_in->set_bit(source);

		Weight upper_bound = init(source, sink); // record current upper_bound
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			Weight tmp = forward(incoming_index_value_source[h_i], outgoing_index_value_sink[h_i]);
			if (merge(tmp, upper_bound) != upper_bound)
			{
				upper_bound = tmp;
			}
		}

		query_engine.template process_vertices<VertexId>(
			[&](VertexId vtx)
			{
				value_source[vtx] = init(source, vtx);
				for (int h_i = 0; h_i < hubs; h_i++)
				{
					Weight relax_val = forward(incoming_index_value_source[h_i], outgoing_index_value[vtx][h_i]);
					if (merge(relax_val, value_source[vtx]) != value_source[vtx])
					{
						value_source[vtx] = relax_val;
					}
				}
				return 0;
			},
			active_all);
		VertexId active_source_vertices = 1;

		for (int i_i = 0; active_source_vertices > 0; i_i++)
		{
			iteration_num++;
			active_num += active_source_vertices;
			active_on_iteration[i_i] += active_source_vertices;
			active_source_out->clear();
			active_source_vertices = query_engine.template process_edges<VertexId, Weight>(
				[&](VertexId src)
				{
					query_engine.emit(src, value_source[src]);
				},
				[&](VertexId src, Weight msg, Adjlist &outgoing_adj)
				{
					VertexId activated = 0;
					for (auto iter : outgoing_adj)
					{
						if (!iter.get_valid(i))
							continue;
						VertexId dst = iter.nbr;
						Weight relax_val = forward(msg, iter.data);
						if (merge(relax_val, value_source[dst]) != value_source[dst])
						{
							if (monotonic_write(&value_source[dst], relax_val))
							{
								if (dst == sink && merge(value_source[dst], upper_bound) != upper_bound)
								{
									monotonic_write(&upper_bound, value_source[dst]); // update shortest path
								}
								active_source_out->set_bit(dst);
								activated += 1;
							}
						}
					}
					return activated;
				},
				[&](VertexId dst, Adjlist &incoming_adj)
				{
					Weight msg = init(source, sink);
					for (auto iter : incoming_adj)
					{
						if (!iter.get_valid(i))
							continue;
						VertexId src = iter.nbr;
						if (active_source_in->get_bit(src))
						{
							Weight relax_val = forward(value_source[src], iter.data);
							if (merge(relax_val, msg) != msg)
							{
								msg = relax_val;
							}
						}
					}
					if (merge(msg, upper_bound) != upper_bound)
						query_engine.emit(dst, msg); // trim by current shortest path
				},
				[&](VertexId dst, Weight msg)
				{
					if (merge(msg, value_source[dst]) != value_source[dst])
					{
						monotonic_write(&value_source[dst], msg);
						if (dst == sink && merge(value_source[dst], upper_bound) != upper_bound)
						{
							monotonic_write(&upper_bound, value_source[dst]); // update shortest path
						}
						active_source_out->set_bit(dst);
						return 1;
					}
					return 0;
				},
				active_source_in);
			std::swap(active_source_in, active_source_out);

			// share upper bound
			if (graph.partition_id == 0)
			{
				for (int p_i = 1; p_i < graph.partitions; p_i++)
				{
					Weight upper_bound_remote;
					MPI_Status status;
					MPI_Recv(&upper_bound_remote, 1, get_mpi_data_type<Weight>(), p_i, QueryMessage, query_comm, &status);
					upper_bound = merge(upper_bound, upper_bound_remote);
				}
			}
			else
			{
				MPI_Send(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, QueryMessage, query_comm);
			}
			MPI_Bcast(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, query_comm);

			active_source_vertices = query_engine.template process_vertices<VertexId>(
				[&](VertexId vtx)
				{
					if (merge(value_source[vtx], upper_bound) == upper_bound)
					{
						active_source_in->clr_bit(vtx);
						return 0;
					}
					return 1;
				},
				active_source_in);
		}
		graph.dealloc_vertex_array_adhoc(value_source);
		delete active_source_in;
		delete active_source_out;
		return upper_bound;
	}

	Weight compute(VertexId source, VertexId sink, uint64_t &active_num, uint64_t &iteration_num, uint64_t *active_on_iteration)
	{
		///////定义变量
		int i = query_snapshot & 1;
		Weight **outgoing_index_value = outgoing_indexes_value[i];
		Weight **incoming_index_value = incoming_indexes_value[i];

		Weight *outgoing_index_value_source = new Weight[hubs];//沿着出边方向找，从中心节点到源顶点的最短距离
		Weight *incoming_index_value_source = new Weight[hubs];//沿着入边方向找，从中心节点到源顶点的最短距离
		Weight *outgoing_index_value_sink = new Weight[hubs];//沿着出边方向找，从中心节点到目的顶点的最短距离
		Weight *incoming_index_value_sink = new Weight[hubs];//沿着入边方向找，从中心节点到目的顶点的最短距离

		////////初始化
		int source_part = graph.get_partition_id(source);
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			if (source_part == graph.partition_id)
			{
				outgoing_index_value_source[h_i] = outgoing_index_value[source][h_i];
				incoming_index_value_source[h_i] = incoming_index_value[source][h_i];
			}
		}
		MPI_Bcast(outgoing_index_value_source, hubs * sizeof(Weight), MPI_CHAR, source_part, query_comm);
		MPI_Bcast(incoming_index_value_source, hubs * sizeof(Weight), MPI_CHAR, source_part, query_comm);

		int sink_part = graph.get_partition_id(sink);
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			if (sink_part == graph.partition_id)
			{
				outgoing_index_value_sink[h_i] = outgoing_index_value[sink][h_i];
				incoming_index_value_sink[h_i] = incoming_index_value[sink][h_i];
			}
		}
		MPI_Bcast(outgoing_index_value_sink, hubs * sizeof(Weight), MPI_CHAR, sink_part, query_comm);
		MPI_Bcast(incoming_index_value_sink, hubs * sizeof(Weight), MPI_CHAR, sink_part, query_comm);

		///////计算上下界
		Weight upper_bound = init(source, sink); // 对于BFS算法，如果起始点和终点不重合，那么上界的初始值是1000000000
#ifdef INDEX_UPPER
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			Weight tmp = forward(incoming_index_value_source[h_i], outgoing_index_value_sink[h_i]);// 对于BFS算法，这句话将两个参数相加
			if (merge(tmp, upper_bound) != upper_bound)
			{
				upper_bound = tmp;
			}
		}
#endif

#ifdef INDEX_LOWER
		/////下界剪枝并，比较巧妙，当下界不满足就直接终止
		Weight lower_bound = init(source, source);//对于BFS算法，下界的初始值是0
		for (int h_i = 0; h_i < hubs; h_i++)
		{
			//这里的下界剪枝，对应于论文中的公式 𝑄(𝑣 ↦ 𝑑) ⪰ 𝑄(ℎ ↦ 𝑑) ⊖ 𝑄(ℎ ↦ 𝑣)、 𝑄(𝑣 ↦ 𝑑) ⪰ 𝑄(𝑣 ↦ ℎ) ⊖ 𝑄(𝑑 ↦ ℎ) 

			Weight tmp;
			tmp = backward(incoming_index_value_sink[h_i], incoming_index_value_source[h_i]);//对于BFS算法，它首先判断参数1是否小于参数2，如果是，则返回（参数2-参数1），否则返回0
			if (merge(tmp, lower_bound) != tmp)
			{
				lower_bound = tmp;
			}
			tmp = backward(outgoing_index_value_source[h_i], outgoing_index_value_sink[h_i]);
			if (merge(tmp, lower_bound) != tmp)
			{
				lower_bound = tmp;
			}
		}
		if (merge(forward(init(source, source), lower_bound), upper_bound) == upper_bound)
		{ // trim by lower_bound
			return upper_bound;//感觉这里写的不好，不能不满足下界就直接返回上界，这样程序就终止了
		}
#endif


		///////定义、初始化待处理顶点集
		Weight *value_source = graph.template alloc_vertex_array_adhoc<Weight>();
		Weight *value_sink = graph.template alloc_vertex_array_adhoc<Weight>();
		VertexSubset *active_source_in = graph.alloc_vertex_subset();
		VertexSubset *active_source_out = graph.alloc_vertex_subset();
		VertexSubset *active_sink_in = graph.alloc_vertex_subset();
		VertexSubset *active_sink_out = graph.alloc_vertex_subset();
		active_source_in->clear();
		active_source_in->set_bit(source);
		active_sink_in->clear();
		active_sink_in->set_bit(sink);

		// template<typename R>
		// R process_vertices(std::function<R(VertexId)> process, Bitmap * active)
		query_engine.template process_vertices<VertexId>(
			[&](VertexId vtx)
			{
				value_source[vtx] = init(source, vtx);
				value_sink[vtx] = init(vtx, sink);
				return 0;
			},
			active_all);
		VertexId active_source_vertices = 1;
		VertexId active_sink_vertices = 1;

#ifdef BFS_OPT
		int intersect = 0;
#endif
		/////////开始计算
#ifdef BIDIRECTIONAL
		for (int i_i = 0; active_source_vertices > 0 || active_sink_vertices > 0; i_i++)
		{
#else
		for (int i_i = 0; active_source_vertices > 0; i_i++)
		{
#endif		//////////通过上下界剪枝，计算从顶点vtx开始需要激活的顶点数量。
			active_source_vertices = query_engine.template process_vertices<VertexId>(
				[&](VertexId vtx)
				{
#ifdef INDEX_UPPER
#ifdef BIDIRECTIONAL
					if (merge(forward(value_source[vtx], value_source[vtx]), upper_bound) == upper_bound)
					{ // trim by upper_bound
#else
					if (merge(value_source[vtx], upper_bound) == upper_bound)
					{ // trim by upper_bound
#endif
						active_source_in->clr_bit(vtx);
						return 0;
					}
#endif
#ifdef INDEX_LOWER
					/////下界剪枝并，比较巧妙，当下界不满足就直接终止
					Weight lower_bound = init(vtx, vtx);
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						Weight tmp;
						tmp = backward(incoming_index_value_sink[h_i], incoming_index_value[vtx][h_i]);
						if (merge(tmp, lower_bound) != tmp)
						{
							lower_bound = tmp;
						}
						tmp = backward(outgoing_index_value[vtx][h_i], outgoing_index_value_sink[h_i]);
						if (merge(tmp, lower_bound) != tmp)
						{
							lower_bound = tmp;
						}
					}
					if (merge(forward(value_source[vtx], lower_bound), upper_bound) == upper_bound)
					{ // trim by lower_bound
						active_source_in->clr_bit(vtx);
						return 0;
					}
#endif
					return 1;
				},
				active_source_in);
			


			///////前向计算
			if (active_source_vertices > 0)
			{ // forward
				iteration_num++;
				active_num += active_source_vertices;
				active_on_iteration[i_i] += active_source_vertices;
				active_source_out->clear();

				// test no buffer vs message buffer (when single node)
				// active_source_vertices = 0;
				// #pragma omp parallel for reduction(+:active_source_vertices)
				// for (VertexId v_i=0;v_i<graph.vertices;v_i++)
				// {
				// 	if (!active_source_in->get_bit(v_i)) continue;
				// 	VertexId local_active_source_vertices = 0;
				// 	uint32_t i = query_snapshot & 1;
				// 	for (auto iter : query_engine.outgoing_adjlist[v_i]) {
				// 		if (!iter.get_valid(i)) continue;
				// 		VertexId dst = iter.nbr;
				// 		Weight relax_val = forward(value_source[v_i], iter.data);
				// 		if (merge(relax_val, value_source[dst]) != value_source[dst]) {
				// 			if (monotonic_write(&value_source[dst], relax_val)) {
				// 				Weight tmp = forward(value_source[dst], value_sink[dst]);
				// 				if (merge(tmp, upper_bound) != upper_bound) { // update upper_bound
				// 					monotonic_write(&upper_bound, tmp);
				// 				}
				// 				active_source_out->set_bit(dst);
				// 				local_active_source_vertices += 1;
				// 			}
				// 		}
				// 	}
				// 	active_source_vertices += local_active_source_vertices;
				// }

				// template<typename R, typename M>
				// R process_edges_sparse(std::function<void(VertexId)> sparse_signal, std::function<R(VertexId, M, Adjlist&)> sparse_slot, Bitmap * active) 
				
				
				active_source_vertices = query_engine.template process_edges_sparse<VertexId, Weight>(
					[&](VertexId src)
					{
						query_engine.emit(src, value_source[src]);
					},
					[&](VertexId src, Weight msg, Adjlist &outgoing_adj)
					{
						VertexId activated = 0;
						uint32_t i = query_snapshot & 1;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							Weight relax_val = forward(msg, iter.data);
							if (merge(relax_val, value_source[dst]) != value_source[dst])
							{
								if (monotonic_write(&value_source[dst], relax_val))
								{
									Weight tmp = forward(value_source[dst], value_sink[dst]);
									if (merge(tmp, upper_bound) != upper_bound)
									{ // update upper_bound
										monotonic_write(&upper_bound, tmp);//单调写入，只有当tmp小于upper_bound时，才会更新upper_bound
#ifdef BFS_OPT
										intersect = 1;//一旦确定路径已经被发现，所有其他的搜索尝试都可以提前终止，从而节省计算资源。
#endif
									}
									active_source_out->set_bit(dst);
									activated += 1;
								}
							}
						}
						return activated;
					},
					active_source_in);

				std::swap(active_source_in, active_source_out);
				MPI_Datatype dt = get_mpi_data_type<Weight>();
				// share upper bound
				if (graph.partition_id == 0)
				{
					for (int p_i = 1; p_i < graph.partitions; p_i++)
					{
						Weight upper_bound_remote;
						MPI_Status status;
						MPI_Recv(&upper_bound_remote, 1, get_mpi_data_type<Weight>(), p_i, QueryMessage, query_comm, &status);
						upper_bound = merge(upper_bound, upper_bound_remote);
					}
				}
				else
				{
					MPI_Send(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, QueryMessage, query_comm);
				}
				MPI_Bcast(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, query_comm);

#ifdef BFS_OPT
				MPI_Allreduce(MPI_IN_PLACE, &intersect, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
				if (intersect)
				{
					break;
				}
#endif
			}

#ifdef BIDIRECTIONAL
			active_sink_vertices = query_engine.template process_vertices<VertexId>(
				[&](VertexId vtx)
				{
#ifdef INDEX_UPPER
					if (merge(forward(value_sink[vtx], value_sink[vtx]), upper_bound) == upper_bound)
					{ // trim by upper_bound
						active_sink_in->clr_bit(vtx);
						return 0;
					}
#endif
#ifdef INDEX_LOWER
					/////下界剪枝并，比较巧妙，当下界不满足就直接终止
					Weight lower_bound = init(vtx, vtx);
					for (int h_i = 0; h_i < hubs; h_i++)
					{
						Weight tmp;
						tmp = backward(incoming_index_value[vtx][h_i], incoming_index_value_source[h_i]);
						if (merge(tmp, lower_bound) != tmp)
						{
							lower_bound = tmp;
						}
						tmp = backward(outgoing_index_value_source[h_i], outgoing_index_value[vtx][h_i]);
						if (merge(tmp, lower_bound) != tmp)
						{
							lower_bound = tmp;
						}
					}
					if (merge(forward(lower_bound, value_sink[vtx]), upper_bound) == upper_bound)
					{ // trim by lower_bound
						active_sink_in->clr_bit(vtx);
						return 0;
					}
#endif
					return 1;
				},
				active_sink_in);

			if (active_sink_vertices > 0)
			{ // backward
				iteration_num++;
				active_num += active_sink_vertices;
				active_on_iteration[i_i] += active_sink_vertices;
				query_engine.transpose();
				active_sink_out->clear();

				// test no buffer vs message buffer (when single node)
				// active_sink_vertices = 0;
				// #pragma omp parallel for reduction(+:active_sink_vertices)
				// for (VertexId v_i=0;v_i<graph.vertices;v_i++)
				// {
				// 	if (!active_sink_in->get_bit(v_i)) continue;
				// 	VertexId local_active_sink_vertices = 0;
				// 	uint32_t i = query_snapshot & 1;
				// 	for (auto iter : query_engine.outgoing_adjlist[v_i]) {
				// 		if (!iter.get_valid(i)) continue;
				// 		VertexId dst = iter.nbr;
				// 		Weight relax_val = forward(value_sink[v_i], iter.data);
				// 		if (merge(relax_val, value_sink[dst]) != value_sink[dst]) {
				// 			if (monotonic_write(&value_sink[dst], relax_val)) {
				// 				Weight tmp = forward(value_source[dst], value_sink[dst]);
				// 				if (merge(tmp, upper_bound) != upper_bound) { // update upper_bound
				// 					monotonic_write(&upper_bound, tmp);
				// 				}
				// 				active_sink_out->set_bit(dst);
				// 				local_active_sink_vertices += 1;
				// 			}
				// 		}
				// 	}
				// 	active_sink_vertices += local_active_sink_vertices;
				// }

				active_sink_vertices = query_engine.template process_edges_sparse<VertexId, Weight>(
					[&](VertexId src)
					{
						query_engine.emit(src, value_sink[src]);
					},
					[&](VertexId src, Weight msg, Adjlist &outgoing_adj)
					{
						VertexId activated = 0;
						uint32_t i = query_snapshot & 1;
						for (auto iter : outgoing_adj)
						{
							if (!iter.get_valid(i))
								continue;
							VertexId dst = iter.nbr;
							Weight relax_val = forward(msg, iter.data);
							if (merge(relax_val, value_sink[dst]) != value_sink[dst])
							{
								if (monotonic_write(&value_sink[dst], relax_val))
								{
									Weight tmp = forward(value_source[dst], value_sink[dst]);
									if (merge(tmp, upper_bound) != upper_bound)
									{ // update upper_bound
										monotonic_write(&upper_bound, tmp);
#ifdef BFS_OPT
										intersect = 1;//一旦确定路径已经被发现，所有其他的搜索尝试都可以提前终止，从而节省计算资源。
#endif
									}
									active_sink_out->set_bit(dst);
									activated += 1;
								}
							}
						}
						return activated;
					},
					active_sink_in);

				query_engine.transpose();
				std::swap(active_sink_in, active_sink_out);
				MPI_Datatype dt = get_mpi_data_type<Weight>();
				// share upper bound
				if (graph.partition_id == 0)
				{
					for (int p_i = 1; p_i < graph.partitions; p_i++)
					{
						Weight upper_bound_remote;
						MPI_Status status;
						MPI_Recv(&upper_bound_remote, 1, get_mpi_data_type<Weight>(), p_i, QueryMessage, query_comm, &status);
						upper_bound = merge(upper_bound, upper_bound_remote);
					}
				}
				else
				{
					MPI_Send(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, QueryMessage, query_comm);
				}
				MPI_Bcast(&upper_bound, 1, get_mpi_data_type<Weight>(), 0, query_comm);

#ifdef BFS_OPT
				MPI_Allreduce(MPI_IN_PLACE, &intersect, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
				if (intersect)
				{
					break;
				}
#endif
			}

#endif
		}

		graph.dealloc_vertex_array_adhoc(value_source);
		graph.dealloc_vertex_array_adhoc(value_sink);
		delete active_source_in;
		delete active_source_out;
		delete active_sink_in;
		delete active_sink_out;
		return upper_bound;
	}
};

#endif
