#include <stdio.h>
#include <stdlib.h>
#include <queue>

#include "core/scheduler.hpp"
#include "core/filesystem.hpp"
using namespace std;

typedef float Weight;
typedef float EdgeData;

class Init {
public:
    Weight operator() (VertexId src, VertexId vid) {
	    return src == vid ? 0 : 1e9;
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
    Scheduler<Weight, EdgeData, Init, Forward, Backward, Merge> scheduler(argc, argv, false);
#ifdef UPDATE
    scheduler.work_update();
#else
    scheduler.work_static();
#endif
    return 0;
}