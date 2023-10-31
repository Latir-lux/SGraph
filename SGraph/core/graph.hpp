#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <malloc.h>
#include <sys/mman.h>
#include <numa.h>
#include <omp.h>

#include <string>
#include <vector>
#include <thread> 
#include <mutex>
#include <functional>
#include <thread>

#include "core/storage.hpp"
#include "core/atomic.hpp"
#include "core/bitmap.hpp"
#include "core/constants.hpp"
#include "core/filesystem.hpp"
#include "core/mpi.hpp"
#include "core/time.hpp"
#include "core/type.hpp"

enum MessageTag
{
	GatherVertexArray,
	QueryMessage,
	IndexMessage
};
template <typename EdgeData>
class Graph
{
public:  
	using Adjlist = std::vector<AdjEdge<EdgeData>>;//AdjEdge是一个结构体，包含了顶点id和边数据,Adjlist是一个vector，里面存放的是AdjEdge

	int partition_id;
	int partitions;

	bool symmetric;
	VertexId vertices;//记录图中顶点的总数量
	EdgeId edge_num;

	uint32_t &query_snapshot;//记录哪个快照处于稳定态,即可以被查询,它的值要么是0,要么是1

	VertexId *partition_offset; // 记录每个分区的偏移量,间接记录分区的顶点数。这是作者的惯用写法，先定义一个指针，然后在构造函数中为其分配内存
	VertexId owned_vertices;//记录本分区的顶点数,借助partition_offset可以计算出来

	VertexId *out_degree; // 记录每个顶点的出度
	VertexId *in_degree;  // 记录每个顶点的入度

	Storage<EdgeData> **outgoing_storage; //这种结构常常用于动态数组或者二维数组的动态分配。在这种情况下，outgoing_storage可能被用作一个二维数组的首地址，而每个outgoing_storage[i]是一个（Storage<EdgeData>*）指针，指向一个Storage<EdgeData>类型的数组。
	Storage<EdgeData> **incoming_storage; // 这种结构常常用于动态数组或者二维数组的动态分配。在这种情况下，incoming_storage可能被用作一个二维数组的首地址，而每个incoming_storage[i]是一个（Storage<EdgeData>*）指针，指向一个Storage<EdgeData>类型的数组。

	std::mutex *outgoing_mutex; // 一个指针, 指向一个互斥锁数组,数组中每个元素都是一个锁
	std::mutex *incoming_mutex; // 一个指针, 指向一个互斥锁数组,数组中每个元素都是一个锁



/*
VertexSubset是对Bitmap的重命名

VertexSubset *: 这意味着我们正在声明一个指向VertexSubset的指针。

outgoing_active[2]: 这意味着outgoing_active是一个包含2个元素的数组，每个元素都是一个VertexSubset *。

在此声明之后，outgoing_active 的布局可能如下：

outgoing_active[0] 是一个指向VertexSubset的指针。
outgoing_active[1] 也是一个指向VertexSubset的指针。
*/
	VertexSubset *outgoing_active[2];//VertexSubset是对Bitmap的重命名
	VertexSubset *incoming_active[2];//VertexSubset是对Bitmap的重命名
	
	Graph(VertexId _vertices, bool _symmetric, uint32_t &_query_snapshot) : query_snapshot(_query_snapshot)
	{
		vertices = _vertices;
		symmetric = _symmetric;
		edge_num = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &partition_id);//获取当前进程的rank
		MPI_Comm_size(MPI_COMM_WORLD, &partitions);//获取进程总数

		partition_offset = new VertexId[partitions + 1];
		partition_offset[0] = 0;
		for (int p_i = 0; p_i < partitions; p_i++)
		{//划分分区,并将分区的偏移量存储在partition_offset中
			partition_offset[p_i + 1] = vertices / partitions * (p_i + 1) / PAGESIZE * PAGESIZE;//先除以PAGESIZE，然后再乘以PAGESIZE，这样就可以得到一个PAGESIZE的倍数,这个数小于等于vertices
			if (p_i == partitions - 1)
			{
				partition_offset[p_i + 1] = vertices;
			}
		}
		owned_vertices = partition_offset[partition_id + 1] - partition_offset[partition_id];//owned_vertices记录本分区的顶点数
		if (symmetric)
		{//如果是对称图，那么 in_degree 和 out_degree 应该是一样的,同理incoming_storage 和 outgoing_storage 也应该是一样的。
			in_degree = out_degree = new VertexId[vertices];
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				out_degree[v_i] = 0;
			}
			incoming_storage = outgoing_storage = new Storage<EdgeData> *[vertices];//new 为一个数组分配堆上的内存,这个数据的每个元素都是一个Storage<EdgeData>类型的指针,数组共有vertices个元素。返回的是数组的首地址。所以outgoing_storage是一个指针，指向一个Storage<EdgeData>类型的数组。incoming_storage同理。
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				outgoing_storage[v_i] = new Storage<EdgeData>(query_snapshot);//根据前面的定义outgoing_storage是一个指针，指向一个Storage<EdgeData>类型的数组，outgoing_storage[v_i]是一个指针类型，指向一个Storage<EdgeData>类型的对象。
			}
			incoming_mutex = outgoing_mutex = new std::mutex[vertices];
			for (int i = 0; i < 2; i++)
			{
				incoming_active[i] = outgoing_active[i] = alloc_vertex_subset();
			}
		}
		else
		{
			out_degree = new VertexId[vertices];
			in_degree = new VertexId[vertices];
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				out_degree[v_i] = 0;
				in_degree[v_i] = 0;
			}
			outgoing_storage = new Storage<EdgeData> *[vertices];
			incoming_storage = new Storage<EdgeData> *[vertices];
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				outgoing_storage[v_i] = new Storage<EdgeData>(query_snapshot);
				incoming_storage[v_i] = new Storage<EdgeData>(query_snapshot);
			}
			outgoing_mutex = new std::mutex[vertices];
			incoming_mutex = new std::mutex[vertices];
			for (int i = 0; i < 2; i++)
			{
				outgoing_active[i] = alloc_vertex_subset();
				incoming_active[i] = alloc_vertex_subset();
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// allocate a vertex array
	template <typename T>
	T *alloc_vertex_array()
	{
		T *array = (T *)malloc(vertices * sizeof(T));
		assert(array != NULL);
		return array;
	}

	// deallocate a vertex array
	template <typename T>
	void dealloc_vertex_array(T *array)
	{
		delete array;
	}

	// allocate a adhoc vertex array (MAKE SURE THAT YOU ONLY VISIT OWNED VERTICES)
	template <typename T>
	T *alloc_vertex_array_adhoc()
	{
		T *array = new T[owned_vertices];
		assert(array != NULL);
		array -= partition_offset[partition_id];//这里将指针位置偏移了partition_offset[partition_id]个位置，是一种巧妙的写法。这样在访问本分区的顶点的时候，就不需要将全局ID转化为局部ID。但是这样做没有边界保护，必须确保分区只访问本分区的顶点，否则可能发生数组越界。
		return array;
	}

	// deallocate a adhoc vertex array
	template <typename T>
	void dealloc_vertex_array_adhoc(T *array)
	{
		array += partition_offset[partition_id];//将array恢复到new返回的原始值
		delete array;
	}

	// allocate a adhoc 2-dimention vertex array (MAKE SURE THAT YOU ONLY VISIT OWNED VERTICES)
	// second dimention has a length of hubs
	template <typename T>
	auto alloc_vertex_array_adhoc_perhub()
	{
		T(*array)
		[hubs] = new T[owned_vertices][hubs];//T(*array)[hubs] 是一个指针声明，其中 array 是一个指向 T[hubs] 类型数组的指针。这意味着 array 指向的是一个包含 hubs 个 T 类型元素的数组。d等式左边，变量名为array，类型为T[hub]（一个数组，数组共hub个元素，每个元素的类型为T)
		assert(array != NULL);
		array -= partition_offset[partition_id];
		return array;
	}

	// deallocate a adhoc 2-dimention vertex array
	// second dimention has a length of hubs
	template <typename T>
	void dealloc_vertex_array_adhoc_perhub(T *array)
	{
		array += partition_offset[partition_id];
		delete[] array;
	}

	// fill a vertex array with a specific value
	template <typename T>
	void fill_vertex_array(T *array, T value)
	{
#pragma omp parallel for
		for (VertexId v_i = partition_offset[partition_id]; v_i < partition_offset[partition_id + 1]; v_i++)
		{
			array[v_i] = value;
		}
	}

	// gather a vertex array
	template <typename T>
	void gather_vertex_array(T *array, int root)//这个函数使用MPI（Message Passing Interface）来收集分布在不同进程上的顶点数组数据。
	{
		if (partition_id != root)//如果当前进程（或称为分区）的ID不是root，那么它应该发送它的数据到root进程。
		{
			MPI_Send(array + partition_offset[partition_id], (uint64_t)sizeof(T) * owned_vertices, MPI_CHAR, root, GatherVertexArray, MPI_COMM_WORLD);
		}
		else
		{
			for (int i = 0; i < partitions; i++)
			{
				if (i == partition_id)
					continue;
				MPI_Status recv_status;
				MPI_Recv(array + partition_offset[i], (uint64_t)sizeof(T) * (partition_offset[i + 1] - partition_offset[i]), MPI_CHAR, i, GatherVertexArray, MPI_COMM_WORLD, &recv_status);//MPI_Recv的参数含义：接收数据的起始地址，接收数据的长度，数据类型，发送方的rank，tag，通信域，状态。
				int length;
				MPI_Get_count(&recv_status, MPI_CHAR, &length);
				assert(length == sizeof(T) * (partition_offset[i + 1] - partition_offset[i]));
			}
		}
	}

	// allocate a vertex subset
	VertexSubset *alloc_vertex_subset()
	{
		return new VertexSubset(vertices);
	}

	int get_partition_id(VertexId v_i)
	{
		for (int i = 0; i < partitions; i++)
		{
			if (v_i >= partition_offset[i] && v_i < partition_offset[i + 1])
			{
				return i;
			}
		}
		assert(false);
	}

	// used to evaluate larger graphs with limited memory
	// init load
	void init_get_degree(EdgeUnit<EdgeData> *edges, long length)
	{//一次性读取一部分边，edges存放边数据，length表示本次读取的边长度
#pragma omp parallel for
		for (long i = 0; i < length; i++)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				write_add(&out_degree[src], (VertexId)1);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				write_add(&in_degree[dst], (VertexId)1);
			}
		}
	}

	// used to evaluate larger graphs with limited memory
	// init load
	void init_load_edges_visible(EdgeUnit<EdgeData> *edges, long length)
	{
		edge_num += length;

#pragma omp parallel for
		for (long i = 0; i < length; i++)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);//使用互斥量来保护存储区域，并避免并发问题。然后创建一个邻接边并将其设置为有效，最后将该邻接边添加到源点的出度存储列表中。
				AdjEdge<EdgeData> adjedge(dst, edges[i].edge_data);
				adjedge.set_valid(0);
				adjedge.set_valid(1);//这两行将valid变量的前两位都设置为1。这意味着这条边在两个不同的snapshot版本（即版本0和版本1）上都被标记为有效。
				outgoing_storage[src]->adjlist.emplace_back(adjedge);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				AdjEdge<EdgeData> adjedge(src, edges[i].edge_data);
				adjedge.set_valid(0);
				adjedge.set_valid(1);
				incoming_storage[dst]->adjlist.emplace_back(adjedge);
			}
		}
	}

	// used to evaluate larger graphs with limited memory
	// init load
	void init_load_edges_invisible(EdgeUnit<EdgeData> *edges, long length)
	{
		edge_num += length;

#pragma omp parallel for
		for (long i = 0; i < length; i++)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				AdjEdge<EdgeData> adjedge(dst, edges[i].edge_data);
				adjedge.clr_valid(0);
				adjedge.clr_valid(1);
				outgoing_storage[src]->adjlist.emplace_back(adjedge);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				AdjEdge<EdgeData> adjedge(src, edges[i].edge_data);
				adjedge.clr_valid(0);
				adjedge.clr_valid(1);
				incoming_storage[dst]->adjlist.emplace_back(adjedge);
			}
		}
	}

	// used to evaluate larger graphs with limited memory
	// update
	void update_edges_visible_to_invisible(EdgeUnit<EdgeData> *edges, long length)
	{
#pragma omp parallel for
		for (long i = 0; i < length; i++)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				AdjEdge<EdgeData> adjedge(dst, edges[i].edge_data);
				for (auto &iter : outgoing_storage[src]->adjlist)
				{
					if (adjedge == iter)
					{
						iter.clr_valid(0);
						iter.clr_valid(1);
						break;
					}
				}
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				AdjEdge<EdgeData> adjedge(src, edges[i].edge_data);
				for (auto &iter : incoming_storage[dst]->adjlist)
				{
					if (adjedge == iter)
					{
						iter.clr_valid(0);
						iter.clr_valid(1);
						break;
					}
				}
			}
		}
	}

	// used to evaluate larger graphs with limited memory
	// update
	void update_edges_invisible_to_visible(EdgeUnit<EdgeData> *edges, long length)
	{
#pragma omp parallel for
		for (long i = 0; i < length; i++)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				AdjEdge<EdgeData> adjedge(dst, edges[i].edge_data);
				for (auto &iter : outgoing_storage[src]->adjlist)
				{
					if (adjedge == iter)
					{
						iter.set_valid(0);
						iter.set_valid(1);
						break;
					}
				}
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				AdjEdge<EdgeData> adjedge(src, edges[i].edge_data);
				for (auto &iter : incoming_storage[dst]->adjlist)
				{
					if (adjedge == iter)
					{
						iter.set_valid(0);
						iter.set_valid(1);
						break;
					}
				}
			}
		}
	}

	void add_edges(std::vector<EdgeUnit<EdgeData>> &edges)
	{
		uint32_t i = query_snapshot & 1;
		edge_num += edges.size();
#pragma omp parallel for
		for (auto iter : edges)
		{
			VertexId src = iter.src;
			VertexId dst = iter.dst;
			write_add(&out_degree[src], (VertexId)1);
			write_add(&in_degree[dst], (VertexId)1);
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, iter.edge_data), true);
				outgoing_active[i]->set_bit(src);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				incoming_storage[dst]->update(AdjEdge<EdgeData>(src, iter.edge_data), true);
				incoming_active[i]->set_bit(dst);
			}
		}
	}

	void del_edges(std::vector<EdgeUnit<EdgeData>> &edges)
	{
		uint32_t i = query_snapshot & 1;
		edge_num -= edges.size();
#pragma omp parallel for
		for (auto iter : edges)
		{
			VertexId src = iter.src;
			VertexId dst = iter.dst;
			write_sub(&out_degree[src], (VertexId)1);
			write_sub(&in_degree[dst], (VertexId)1);
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, iter.edge_data), false);
				outgoing_active[i]->set_bit(src);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				incoming_storage[dst]->update(AdjEdge<EdgeData>(src, iter.edge_data), false);
				incoming_active[i]->set_bit(dst);
			}
		}
	}

	void add_edge_atomic(EdgeUnit<EdgeData> &edge)
	{
		uint32_t i = query_snapshot & 1;
		edge_num += 1;
		VertexId src = edge.src;
		VertexId dst = edge.dst;
		write_add(&out_degree[src], (VertexId)1);
		write_add(&in_degree[dst], (VertexId)1);
		if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
		{
			std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
			outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, edge.edge_data), true);
			outgoing_active[i]->set_bit(src);
		}
		if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
		{
			std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
			incoming_storage[dst]->update(AdjEdge<EdgeData>(src, edge.edge_data), true);
			incoming_active[i]->set_bit(dst);
		}
	}

	void del_edge_atomic(EdgeUnit<EdgeData> &edge)
	{
		uint32_t i = query_snapshot & 1;
		edge_num -= 1;
		VertexId src = edge.src;
		VertexId dst = edge.dst;
		write_sub(&out_degree[src], (VertexId)1);
		write_sub(&in_degree[dst], (VertexId)1);
		if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
		{
			std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
			outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, edge.edge_data), false);
			outgoing_active[i]->set_bit(src);
		}
		if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
		{
			std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
			incoming_storage[dst]->update(AdjEdge<EdgeData>(src, edge.edge_data), false);
			incoming_active[i]->set_bit(dst);
		}
	}

	void step()
	{
		uint32_t i = query_snapshot & 1;
#pragma omp parallel for
		for (VertexId v_i = 0; v_i < vertices; v_i++)
		{
			if (outgoing_active[0]->get_bit(v_i) || outgoing_active[1]->get_bit(v_i))
			{
				outgoing_storage[v_i]->step();
			}
		}
		outgoing_active[i ^ 1]->clear();
		if (!symmetric)
		{
#pragma omp parallel for
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				if (incoming_active[0]->get_bit(v_i) || incoming_active[1]->get_bit(v_i))
				{
					incoming_storage[v_i]->step();
				}
			}
			incoming_active[i ^ 1]->clear();
		}
	}
};

#endif
