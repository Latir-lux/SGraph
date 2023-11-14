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
	using Adjlist = std::vector<AdjEdge<EdgeData>>;

	int partition_id;
	int partitions;

	bool symmetric;
	VertexId vertices;
	EdgeId edge_num;

	uint32_t &query_snapshot;

	VertexId *partition_offset; // VertexId [partitions+1]
	VertexId owned_vertices;

	VertexId *out_degree; // VertexId [vertices];
	VertexId *in_degree;  // VertexId [vertices];

	Storage<EdgeData> **outgoing_storage; // Storage * [vertices];
	Storage<EdgeData> **incoming_storage; // Storage * [vertices];

	std::mutex *outgoing_mutex; // mutex [vertices]
	std::mutex *incoming_mutex; // mutex [vertices]

	VertexSubset *outgoing_active[2];
	VertexSubset *incoming_active[2];

	Graph(VertexId _vertices, bool _symmetric, uint32_t &_query_snapshot) : query_snapshot(_query_snapshot)
	{
		vertices = _vertices;
		symmetric = _symmetric;
		edge_num = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &partition_id);
		MPI_Comm_size(MPI_COMM_WORLD, &partitions);

		partition_offset = new VertexId[partitions + 1];
		partition_offset[0] = 0;
		for (int p_i = 0; p_i < partitions; p_i++)
		{
			partition_offset[p_i + 1] = vertices / partitions * (p_i + 1) / PAGESIZE * PAGESIZE;
			if (p_i == partitions - 1)
			{
				partition_offset[p_i + 1] = vertices;
			}
		}
		owned_vertices = partition_offset[partition_id + 1] - partition_offset[partition_id];
		if (symmetric)
		{
			in_degree = out_degree = new VertexId[vertices];
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				out_degree[v_i] = 0;
			}
			incoming_storage = outgoing_storage = new Storage<EdgeData> *[vertices];
			for (VertexId v_i = 0; v_i < vertices; v_i++)
			{
				outgoing_storage[v_i] = new Storage<EdgeData>(query_snapshot);
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
		array -= partition_offset[partition_id];
		return array;
	}

	// deallocate a adhoc vertex array
	template <typename T>
	void dealloc_vertex_array_adhoc(T *array)
	{
		array += partition_offset[partition_id];
		delete array;
	}

	// allocate a adhoc 2-dimention vertex array (MAKE SURE THAT YOU ONLY VISIT OWNED VERTICES)
	// second dimention has a length of hubs
	template <typename T>
	auto alloc_vertex_array_adhoc_perhub()
	{
		T(*array)
		[hubs] = new T[owned_vertices][hubs];
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
	void gather_vertex_array(T *array, int root)
	{
		if (partition_id != root)
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
				MPI_Recv(array + partition_offset[i], (uint64_t)sizeof(T) * (partition_offset[i + 1] - partition_offset[i]), MPI_CHAR, i, GatherVertexArray, MPI_COMM_WORLD, &recv_status);
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
	{
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
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				AdjEdge<EdgeData> adjedge(dst, edges[i].edge_data);
				adjedge.set_valid(0);
				adjedge.set_valid(1);
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
		for(int i = 0; i < edges.size(); ++i)
		//for (auto iter : edges)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			write_add(&out_degree[src], (VertexId)1);
			write_add(&in_degree[dst], (VertexId)1);
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, edges[i].edge_data), true);
				outgoing_active[i]->set_bit(src);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				incoming_storage[dst]->update(AdjEdge<EdgeData>(src, edges[i].edge_data), true);
				incoming_active[i]->set_bit(dst);
			}
		}
	}

	void del_edges(std::vector<EdgeUnit<EdgeData>> &edges)
	{
		uint32_t i = query_snapshot & 1;
		edge_num -= edges.size();
#pragma omp parallel for
		for(int i = 0; i < edges.size(); ++i)
		//for (auto iter : edges)
		{
			VertexId src = edges[i].src;
			VertexId dst = edges[i].dst;
			write_sub(&out_degree[src], (VertexId)1);
			write_sub(&in_degree[dst], (VertexId)1);
			if (dst >= partition_offset[partition_id] && dst < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(outgoing_mutex[src]);
				outgoing_storage[src]->update(AdjEdge<EdgeData>(dst, edges[i].edge_data), false);
				outgoing_active[i]->set_bit(src);
			}
			if (src >= partition_offset[partition_id] && src < partition_offset[partition_id + 1])
			{
				std::lock_guard<std::mutex> lock(incoming_mutex[dst]);
				incoming_storage[dst]->update(AdjEdge<EdgeData>(src, edges[i].edge_data), false);
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
