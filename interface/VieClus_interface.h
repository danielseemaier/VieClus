#ifndef VIECLUS_VIECLUS_H
#define VIECLUS_VIECLUS_H

namespace VieClus {

struct Graph {
	int n;
	int *xadj;
	int *adjncy;
	int *vwgt;
	int *adjwgt;
	int *clustering; // Nullable
};

__attribute__((visibility("default"))) void setup(int *argc, char ***argv);

__attribute__((visibility("default"))) double run_default(Graph graph, int time_limit, int seed, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) double run_shallow(Graph graph, int time_limit, int seed, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) double run_shallow_no_lp(Graph graph, int time_limit, int seed, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) void teardown();

}

#endif // VIECLUS_VIECLUS_H
