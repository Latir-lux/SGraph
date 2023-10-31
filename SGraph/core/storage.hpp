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
	{ // preprocess is included  update函数的核心思想是延迟对adjlist的实际修改，只是记录要进行的操作。实际的修改将在step函数中完成。

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
		{//使用线性搜索或哈希映射尝试找到给定的边。
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
			{//如果边未找到并且操作是添加，那么在adjlist的末尾添加该边，并在adjlist_map中（如果存在的话）为其设置一个新的映射。
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

		buffer[i].emplace_back(buffer_unit(edge, add, pos));//在当前步骤的buffer中为这个边添加一个新的buffer_unit，其中包含边、操作类型（添加/删除）以及边在adjlist中的位置
	}

	void step()
	{//step函数的目的是应用在buffer中记录的所有更改到adjlist。
		uint32_t i = query_snapshot & 1;//通过切换query_snapshot的最低有效位来确定哪个buffer包含上一步的更改，哪个buffer用于当前步骤的更改。
		for (auto &iter : buffer[i ^ 1])
		{
			if (iter.add == true)//遍历上一步的buffer，并对adjlist中的边应用更改：如果buffer_unit中的add为true，则设置该边为有效；如果为false，则设置该边为无效
			{
				adjlist[iter.pos].set_valid(i);
			}
			else
			{
				adjlist[iter.pos].clr_valid(i);
			}
		}
		buffer[i ^ 1].clear();//清空上一步的buffer，以为新的更改释放空间。
		for (auto &iter : buffer[i])
		{//遍历当前步骤的buffer，并对adjlist中的边应用相同的更改
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
