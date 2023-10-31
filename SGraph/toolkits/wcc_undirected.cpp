#include <stdio.h>
#include <stdlib.h>
#include <queue>

#include "core/scheduler.hpp"
#include "core/filesystem.hpp"
using namespace std;

typedef bool Weight;
typedef Empty<bool> EdgeData;
template <> bool AdjEdge <EdgeData>::data = true;
template <> bool EdgeUnit <EdgeData>::edge_data = true;

class Init {
public:
    Weight operator() (VertexId src, VertexId vid) {
	    return src == vid ? true : false;
    }
};
class Forward {
public:
    Weight operator() (Weight ab, Weight bc) {
	    return ab & bc;
    }
};
class Backward {
public:
    Weight operator() (Weight ab, Weight ac) {
	    return (ab && !ac) ? false : true;
    }
};
class Merge {
public:
    Weight operator() (Weight a, Weight b) {
	    return a | b;
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