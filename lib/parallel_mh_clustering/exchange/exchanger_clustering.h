/******************************************************************************
 * exchanger_clustering.h
 *
 * Source of VieClus -- Vienna Graph Clustering 
 ******************************************************************************
 * Copyright (C) 2017 Sonja Biedermann, Christian Schulz and Bernhard Schuster
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef EXCHANGER_YPB6QKNL
#define EXCHANGER_YPB6QKNL

#include <mpi.h>

#include "data_structure/graph_access.h"
#include "parallel_mh_clustering/population_clustering.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"

class exchanger_clustering {
public:
        exchanger_clustering( MPI_Comm communicator );
        virtual ~exchanger_clustering();

        void diversify_population_clustering( PartitionConfig & config, graph_access & G, population_clustering & island, bool replace );
        void quick_start( PartitionConfig & config,  graph_access & G, population_clustering & island );
        void push_best( PartitionConfig & config,  graph_access & G, population_clustering & island );
        void recv_incoming( PartitionConfig & config,  graph_access & G, population_clustering & island );

private:
        void exchange_individum(const PartitionConfig & config, 
                                graph_access & G, 
                                int & from, 
                                int & rank, 
                                int & to, 
                                Individuum & in, Individuum & out);

        std::vector< int* >          m_partition_map_buffers;
        std::vector< MPI_Request* > m_request_pointers;
        std::vector<bool>            m_allready_send_to;

        double m_prev_best_objective;
        int m_max_num_pushes;
        int m_cur_num_pushes;

        MPI_Comm m_communicator;

        quality_metrics m_qm;
};



#endif /* end of include guard: EXCHANGER_YPB6QKNL */
