#ifndef FILESYSTEM_HPP
#define FILESYSTEM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <mpi.h>

#include "core/type.hpp"
#include "core/constants.hpp"

inline long file_size(std::string filename)
{
	struct stat st;
	assert(stat(filename.c_str(), &st) == 0);
	return st.st_size;
}

inline long file_size(std::string filename, int partition_id)
{
	long ret;
	if (partition_id == 0)
	{
		struct stat st;
		assert(stat(filename.c_str(), &st) == 0);
		ret = st.st_size;
	}
	MPI_Bcast(&ret, 2, MPI_INT, 0, MPI_COMM_WORLD);
	return ret;
}

template <typename EdgeData>
void get_edge_vector(std::string path, std::vector<EdgeUnit<EdgeData>> &edges)
{
	edges.clear();
	int fin = open(path.c_str(), O_RDWR);
	long total_bytes = file_size(path.c_str());
	size_t edge_unit_size = sizeof(EdgeUnit<EdgeData>);
	long total_edges = total_bytes / edge_unit_size;
	long read_edges = 0;
	EdgeUnit<EdgeData> *buffer = new EdgeUnit<EdgeData>[CHUNKSIZE];
	while (read_edges < total_edges)
	{
		long current_read_edges;
		if (total_edges - read_edges > CHUNKSIZE)
		{
			current_read_edges = read(fin, buffer, edge_unit_size * CHUNKSIZE) / edge_unit_size;
		}
		else
		{
			current_read_edges = read(fin, buffer, edge_unit_size * (total_edges - read_edges)) / edge_unit_size;
		}
		for (int i = 0; i < current_read_edges; i++)
		{
			edges.emplace_back(buffer[i]);
		}
		read_edges += current_read_edges;
	}
	close(fin);
}

template <typename EdgeData>
void get_edge_vector(std::string path, EdgeUnit<EdgeData> *edges, long offset_begin, long offset_end, int partition_id)
{
	long total_edges = offset_end - offset_begin;
	size_t edge_unit_size = sizeof(EdgeUnit<EdgeData>);
	if (partition_id == 0)
	{
		int fin = open(path.c_str(), O_RDWR);
		assert(fin != -1);
		lseek(fin, offset_begin * edge_unit_size, SEEK_SET);
		long read_edges = 0;
		while (read_edges < total_edges)
		{
			long current_read_edges;
			if (total_edges - read_edges > CHUNKSIZE)
			{
				current_read_edges = read(fin, edges + read_edges, edge_unit_size * CHUNKSIZE) / edge_unit_size;
			}
			else
			{
				current_read_edges = read(fin, edges + read_edges, edge_unit_size * (total_edges - read_edges)) / edge_unit_size;
			}
			read_edges += current_read_edges;
		}
		close(fin);
	}
	MPI_Bcast(edges, total_edges * edge_unit_size, MPI_CHAR, 0, MPI_COMM_WORLD);
}

#endif
