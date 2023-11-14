#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <cstdio>
#include <mutex>
#include <sparsehash/dense_hash_map>
#include <vector>
#include "core/atomic.hpp"
#include "core/type.hpp"
#include "core/time.hpp"

// used to evaluate larger graphs with limited memory
#ifdef STATIC
template <typename EdgeData>
class Storage
{
public:
	using Adjlist = std::vector<AdjEdge<EdgeData>>;

	Adjlist adjlist;

	Storage(uint32_t &_query_snapshot)
	{
	}

	void update(AdjEdge<EdgeData> edge, bool add)
	{ // preprocess is included
	}

	void step()
	{
	}
};

#else

template <typename EdgeData>
class Storage
{
public:
	using Adjlist = std::vector<AdjEdge<EdgeData>>;

	Adjlist adjlist;
	google::dense_hash_map<AdjEdge<EdgeData>, EdgeId> *adjlist_map;
	bool has_map;
	int adjlist_size;
	uint32_t &query_snapshot;

	struct buffer_unit
	{
		AdjEdge<EdgeData> edge;
		bool add;
		int pos;
		buffer_unit(AdjEdge<EdgeData> _edge, bool _add, int _pos)
		{
			edge = _edge;
			add = _add;
			pos = _pos;
		}
	};
	std::vector<buffer_unit> buffer[2];

	Storage(uint32_t &_query_snapshot) : query_snapshot(_query_snapshot)
	{
		has_map = false;
		adjlist_size = 0;
		query_snapshot = 0;
	}

	void update(AdjEdge<EdgeData> edge, bool add)
	{ // preprocess is included

		if (!has_map && adjlist_size >= 512)
		{
			adjlist_map = new google::dense_hash_map<AdjEdge<EdgeData>, EdgeId>;
			AdjEdge<EdgeData> empty_key;
			adjlist_map->set_empty_key(empty_key);
			for (int i = 0; i < adjlist_size; i++)
			{
				adjlist_map->insert(std::make_pair(adjlist[i], i));
			}
			has_map = true;
		}

		uint32_t i = query_snapshot & 1;
		int pos = -1;
		if (!has_map)
		{
			for (int k = 0; k < adjlist_size; k++)
			{
				if (adjlist[k] == edge)
				{
					pos = k;
					break;
				}
			}
		}
		else
		{
			auto iter1 = adjlist_map->find(edge);
			if (iter1 != adjlist_map->end())
			{
				pos = (*iter1).second;
			}
		}
		if (pos == -1)
		{
			if (add == true)
			{
				pos = adjlist_size;
				adjlist_size++;
				adjlist.emplace_back(edge);
				if (has_map)
				{
					adjlist_map->insert(std::make_pair(edge, pos));
				}
			}
			else
			{ // in case edge is added in the same step
				for (auto iter1 : buffer[i])
				{
					if (iter1.edge == edge && iter1.add == true)
					{
						pos = iter1.pos;
						break;
					}
				}
			}
		}

		buffer[i].emplace_back(buffer_unit(edge, add, pos));
	}

	void step()
	{
		uint32_t i = query_snapshot & 1;
		for (auto &iter : buffer[i ^ 1])
		{
			if (iter.add == true)
			{
				adjlist[iter.pos].set_valid(i);
			}
			else
			{
				adjlist[iter.pos].clr_valid(i);
			}
		}
		buffer[i ^ 1].clear();
		for (auto &iter : buffer[i])
		{
			if (iter.add == true)
			{
				if (iter.pos >= adjlist_size)
				{
					// adjlist.emplace_back(iter.edge);
				}
				adjlist[iter.pos].set_valid(i);
			}
			else
			{
				adjlist[iter.pos].clr_valid(i);
			}
		}
	}
};

#endif

#endif
