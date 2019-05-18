#include "VieClus_interface.h"

#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "algorithms/cycle_search.h"
#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

namespace VieClus {

static double _run(const Graph &graph, PartitionConfig &partition_config, int *out_k, int *out_partition_map);

__attribute__((visibility("default"))) void setup(int *argc, char ***argv) {
	MPI_Init(argc, argv);
	omp_set_num_threads(1);
}

__attribute__((visibility("default"))) double run_default(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map) {
	PartitionConfig partition_config;
	configuration cfg;
	cfg.standard(partition_config);
	cfg.strong(partition_config);
	partition_config.seed = seed;
	partition_config.time_limit = time_limit;
	partition_config.k = 1;
	partition_config.bcc_max_lv_iterations = max_lv_iterations;
	partition_config.bcc_upper_cluster_weight = 0;
	return _run(graph, partition_config, out_k, out_partition_map);
}

__attribute__((visibility("default"))) double run_shallow(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map) {
	PartitionConfig partition_config;
	configuration cfg;
	cfg.standard(partition_config);
	cfg.strong(partition_config);
	partition_config.seed = seed;
	partition_config.time_limit = time_limit;
	partition_config.k = 1;
	partition_config.bcc_shallow_coarsening = true;
  partition_config.bcc_max_lv_iterations = max_lv_iterations;
  partition_config.bcc_upper_cluster_weight = 0;
	return _run(graph, partition_config, out_k, out_partition_map);
}

__attribute__((visibility("default"))) double run_shallow_no_lp(Graph graph, int time_limit, int max_lv_iterations, int seed, int *out_k, int *out_partition_map) {
	PartitionConfig partition_config;
	configuration cfg;
	cfg.standard(partition_config);
	cfg.strong(partition_config);
	partition_config.seed = seed;
	partition_config.time_limit = time_limit;
	partition_config.k = 1;
	partition_config.bcc_shallow_coarsening = true;
	partition_config.bcc_no_lp = true;
  partition_config.bcc_max_lv_iterations = max_lv_iterations;
  partition_config.bcc_upper_cluster_weight = 0;
	return _run(graph, partition_config, out_k, out_partition_map);
}

__attribute__((visibility("default"))) double run_shallow_no_lp_size_constrained(Graph graph, int time_limit, int max_lv_iterations, unsigned int upper_cluster_weight, int seed, int *out_k, int *out_partition_map) {
  PartitionConfig partition_config;
  configuration cfg;
  cfg.standard(partition_config);
  cfg.strong(partition_config);
  partition_config.seed = seed;
  partition_config.time_limit = time_limit;
  partition_config.k = 1;
  partition_config.bcc_shallow_coarsening = true;
  partition_config.bcc_no_lp = true;
  partition_config.bcc_max_lv_iterations = max_lv_iterations;
  partition_config.bcc_upper_cluster_weight = upper_cluster_weight;
  return _run(graph, partition_config, out_k, out_partition_map);
}

__attribute__((visibility("default"))) void teardown() {
	MPI_Finalize();
}

static double _run(const Graph &graph, PartitionConfig &partition_config, int *out_k, int *out_partition_map) {
    // prepare graph from Graph data structure
	graph_access G;
	if (graph.vwgt == nullptr || graph.adjwgt == nullptr) {
		G.build_from_metis(graph.n, graph.xadj, graph.adjncy);
	} else {
		G.build_from_metis_weighted(graph.n, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt);
	}

	double initial_modularity = 0.0;
    if (graph.clustering != nullptr) {
        for (NodeID v = 0; v < G.number_of_nodes(); ++v) {
            G.setPartitionIndex(v, graph.clustering[v]);
        }
        G.set_partition_count(G.get_partition_count_compute());
        partition_config.bcc_start_w_singletons = false;

        initial_modularity = ModularityMetric::computeModularity(G);
    }

    // run VieClus clustering algorithm
	parallel_mh_async_clustering mh;
	mh.perform_partitioning(partition_config, G);

	// write results to output parameters
	G.set_partition_count(G.get_partition_count_compute());
	if (out_k != nullptr) *out_k = G.get_partition_count();
	if (out_partition_map != nullptr) {
		for (NodeID v = 0; v < G.number_of_nodes(); ++v) {
			out_partition_map[v] = G.getPartitionIndex(v);
		}
	}

	double final_modularity = ModularityMetric::computeModularity(G);
	if (graph.clustering == nullptr) {
	    std::cout << "[BCCInfo] Started with a singleton clustering" << std::endl;
	    std::cout << "[BCC] final_modularity=" << final_modularity << std::endl;
	} else {
	    std::cout << "[BCCInfo] Improved a given clustering" << std::endl;
	    std::cout << "[BCC] initial_modularity=" << initial_modularity << "; final_modularity=" << final_modularity << std::endl;
	}
	return final_modularity;
}

}