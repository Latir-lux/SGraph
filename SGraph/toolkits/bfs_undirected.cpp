#include <stdio.h>
#include <stdlib.h>
#include <queue>

#include "core/scheduler.hpp"
#include "core/filesystem.hpp"
using namespace std;

typedef int Weight;
typedef Empty<int> EdgeData;
template <> int AdjEdge <EdgeData>::data = 1;
template <> int EdgeUnit <EdgeData>::edge_data = 1;

class Init {
public:
    Weight operator() (VertexId src, VertexId vid) {
	    return src == vid ? 0 : 1000000000;
    }
};
class Forward {
public:
    Weight operator() (Weight ab, Weight bc) {
	    return ab + bc;
    }
};
class Backward {
public:
    Weight operator() (Weight ab, Weight ac) {
	    return ab < ac ? ac - ab : 0;
    }
};
class Merge {
public:
    Weight operator() (Weight a, Weight b) {
	    return a < b ? a : b;
    }
};

int main(int argc, char ** argv){
    Scheduler<Weight, EdgeData, Init, Forward, Backward, Merge> scheduler(argc, argv, true);
#ifdef UPDATE
    scheduler.work_update();
#else
    scheduler.work_static();
#endif
    return 0;
}