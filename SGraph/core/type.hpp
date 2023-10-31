#ifndef TYPE_HPP

#define TYPE_HPP

#include <stdint.h>

template <typename EdgeData>
struct Empty//这是一个空模板结构体,有一个类型参数 EdgeData
{
};

typedef uint32_t VertexId;
const int VertexIdBits = 30;//int类型为32位，VertexIdBits为30，所以VertexId的范围为0~2^30-1
typedef uint64_t EdgeId;//边的类型是uint64_t,点的类型是uint32_t,可能是因为边的数量大于点的数量,所以要采用更大的类型来存储

template <typename EdgeData>
struct EdgeUnit//这是一个模板结构体，有一个类型参数 EdgeData,内部属性:边的起点、终点和边的数据
{
	VertexId src;
	VertexId dst;
	EdgeData edge_data;
};

template <typename EdgeData>
struct EdgeUnit<Empty<EdgeData>>
{//这是 EdgeUnit 的一个特化版本，适用于 Empty<EdgeData> 类型的参数。
	VertexId src;
	VertexId dst;
	static EdgeData edge_data;//edge_data 是静态的：这意味着所有 EdgeUnit<Empty<EdgeData>> 实例都将共享相同的 edge_data 成员。
	EdgeUnit()
	{
		src = -1;
		dst = -1;
	}
	EdgeUnit(EdgeUnit<EdgeData> &a)
	{
		src = a.src;
		dst = a.dst;
	}
};

template <typename EdgeData>
struct AdjEdge
{//这是一个模板结构体，有一个类型参数 EdgeData,内部属性:邻居数,有效性,边的数据
	VertexId nbr : VertexIdBits;
	VertexId valid : 2;
	EdgeData data;
	AdjEdge()
	{
		nbr = (1 << VertexIdBits) - 1;
		memset(&data, 0, sizeof(EdgeData));
		valid = 0;
	}
	AdjEdge(VertexId _nbr, EdgeData _data)
	{
		nbr = _nbr;
		data = _data;
		valid = 0;
	}
	inline bool get_valid(uint32_t snapshot)//valid是一个2位的变量，这里的snapshot是一个0或1的值,get_valid函数返回在snapshot指示的位上的值
	{
		return (valid >> snapshot) & 1;
	}
	inline void set_valid(uint32_t snapshot)//valid是一个2位的变量，这里的snapshot是一个0或1的值，set_valid函数将在snapshot指示的位上的值置为1
	{
		valid |= (1 << snapshot);
	}
	inline void clr_valid(uint32_t snapshot)////valid是一个2位的变量，这里的snapshot是一个0或1的值，clr_valid函数将在snapshot指示的位上的值置为0,其它位不变
	{
		valid &= (1 << (1 ^ snapshot));
	}
	bool operator<(const AdjEdge<EdgeData> &a) const
	{
		if (nbr != a.nbr)
			return nbr < a.nbr;
		return data < a.data;
	}
	bool operator==(const AdjEdge<EdgeData> &a) const
	{
		return (nbr == a.nbr && data == a.data);
	}
};


// 完全特化：所有的模板参数都被具体指定。例如，对于一个模板template <typename T> struct Example;，完全特化可能看起来像template <> struct Example<int> {...};，其中T被明确指定为int。
// 通用模板：模板参数是完全开放的，没有任何约束。使用上面的Example模板，它的通用定义可能是template <typename T> struct Example {...};
// 偏特化：部分模板参数具有某种形式或约束，但仍然保留了某种不确定性。使用上面的Example模板，一个偏特化可能看起来像template <typename T> struct Example<SomeType<T>> {...};，这里的不确定性在于SomeType<T>可以是SomeType<int>, SomeType<double>, SomeType<std::string>等等。
template <typename EdgeData>
struct AdjEdge<Empty<EdgeData>>
{//这是 AdjEdge 的一个特化版本，适用于 Empty<EdgeData> 类型的参数。
	VertexId nbr : VertexIdBits;//这种变量后加冒号的语法是C++中的“位字段”（Bit-field）定义。位字段允许你为整数类型的结构体或联合体成员指定一个确切的宽度（以位为单位）。
	VertexId valid : 2;
	static EdgeData data;
	AdjEdge()
	{
		nbr = (1 << VertexIdBits) - 1;
		valid = 0;
	}
	AdjEdge(VertexId _nbr, EdgeData _data)
	{
		nbr = _nbr;
		valid = 0;
	}
	inline bool get_valid(uint32_t snapshot)
	{
		return (valid >> snapshot) & 1;
	}
	inline void set_valid(uint32_t snapshot)
	{
		valid |= (1 << snapshot);
	}
	inline void clr_valid(uint32_t snapshot)
	{
		valid &= (1 << (1 ^ snapshot));
	}
	bool operator<(const AdjEdge<Empty<EdgeData>> &a) const
	{
		return nbr < a.nbr;
	}
	bool operator==(const AdjEdge<Empty<EdgeData>> &a) const
	{
		return nbr == a.nbr;
	}
};

namespace std
{
	template <typename T>
	struct hash<AdjEdge<T>>
	{//利用邻居数和边数据产生哈希值
		uint64_t operator()(const AdjEdge<T> &e) const
		{
			return (17lu * std::hash<VertexId>()(e.nbr) + std::hash<T>()(e.data));//将 e.nbr 的哈希值乘以一个常数（在这里是 17lu）。乘以一个质数（例如17）是一种常用的哈希组合技术，可以帮助降低不同输入产生相同哈希值的机会。
		}
	};
	template <typename EdgeData>
	struct hash<AdjEdge<Empty<EdgeData>>>
	{//利用邻居数产生哈希值
		uint64_t operator()(const AdjEdge<Empty<EdgeData>> &e) const
		{
			return std::hash<VertexId>()(e.nbr);
		}
	};
}

#endif
