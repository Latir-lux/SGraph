#ifndef TYPE_HPP
#define TYPE_HPP

#include <stdint.h>

template <typename EdgeData>
struct Empty
{
};

typedef uint32_t VertexId;
const int VertexIdBits = 30;
typedef uint64_t EdgeId;

template <typename EdgeData>
struct EdgeUnit
{
	VertexId src;
	VertexId dst;
	EdgeData edge_data;
};

template <typename EdgeData>
struct EdgeUnit<Empty<EdgeData>>
{
	VertexId src;
	VertexId dst;
	static EdgeData edge_data;
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
{
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

template <typename EdgeData>
struct AdjEdge<Empty<EdgeData>>
{
	VertexId nbr : VertexIdBits;
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
	{
		uint64_t operator()(const AdjEdge<T> &e) const
		{
			return (17lu * std::hash<VertexId>()(e.nbr) + std::hash<T>()(e.data));
		}
	};
	template <typename EdgeData>
	struct hash<AdjEdge<Empty<EdgeData>>>
	{
		uint64_t operator()(const AdjEdge<Empty<EdgeData>> &e) const
		{
			return std::hash<VertexId>()(e.nbr);
		}
	};
}

#endif
