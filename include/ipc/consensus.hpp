#pragma once

#include "ipc/consensus_utils.hpp"

// Using full deterministic search
template <class T, class EDGE, class VERTEX>
void simulating_incremental_data(const Config& cfg,
                                 g2o::SparseOptimizer& open_loop_problem,
                                 const std::vector<EDGE*>& loops);

template <class EDGE, class VERTEX> class IPC
{
public :

    IPC(g2o::SparseOptimizer& open_loop_problem, const Config& cfg);
    ~IPC();

    bool agreementCheck(EDGE* loop_candidate);

    const std::vector<std::pair<int, int>>& getMaxConsensusSet() const { return _max_consensus_set; }
    const std::vector<std::pair<int, int>>& getClusters() const { return _clusters; }

private :


    double iou(const std::pair<int, int>& lp1, const std::pair<int, int>& lp2);

    // Return the voters with highest percentage of overlap with the candidate
    // Plus the beginning and end of the cluster
    std::pair<int, int> getTopKVoters(const std::pair<int, int>& candidate, 
                                      const g2o::OptimizableGraph::EdgeSet& eset_loops,
                                      int k, g2o::OptimizableGraph::EdgeSet& voters);

    g2o::SparseOptimizer* _problem;

    std::vector<std::pair<int, int>> _max_consensus_set;
    std::vector<std::pair<int, int>> _clusters;
    std::vector<EDGE*> _odom_edges;
    std::vector<g2o::OptimizableGraph::EdgeSet> _edgesxcl;
    std::vector<g2o::OptimizableGraph::EdgeSet> _edgesxcl_only_loops;

    double _fast_reject_th;
    int _fast_reject_iter_base;
    double _slow_reject_th;
    int _slow_reject_iter_base;
    bool _use_best_k_buddies;
    int _k_buddies;


};