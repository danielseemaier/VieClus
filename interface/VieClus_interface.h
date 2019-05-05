#ifndef VIECLUS_VIECLUS_H
#define VIECLUS_VIECLUS_H

namespace VieClus {

struct Graph {
	int n; //! Number of nodes in the graph
	int *xadj; //! Nodes array
	int *adjncy; //! Edges array
	int *vwgt; //! Node weights
	int *adjwgt; //! Edge weights
	int *clustering; //! Preexisting clustering to be used as initial clustering (nullable)
};

//! Call before using any of the other methods, initializes MPI
__attribute__((visibility("default"))) void setup(int *argc, char ***argv);

/*
 * run_default: run VieClus with default parameters
 * run_shallow: run single-level VieClus, i.e. no contractions etc
 * run_shallow_no_lp: same as run_shallow, but do not shrink the graph using label propagation before running LV
 *
 * Arguments (same for all interface methods):
 * - graph: graph
 * - time_limit: time limit for the evolutionary algorithm, set to 0 to disable the evolutionary algorithm
 * - max_lv_iterations: run at most this many LV iterations, set to 0 to use the normal convergence criteria
 * - seed: seed for random operations
 * - out_k: memory address to store the number of clusters found (nullable)
 * - out_partition_map: memory address to store the clustering (nullable)
 *
 * Return value: modularity of the clustering
 */

__attribute__((visibility("default"))) double run_default(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) double run_shallow(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) double run_shallow_no_lp(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map);

//! Must be the last call to VieClus, finalizes MPI
__attribute__((visibility("default"))) void teardown();

}

#endif // VIECLUS_VIECLUS_H
