#pragma once

#include "ipc/consensus_utils.hpp"

// Using full deterministic search 
void simulating_incremental_data(const Config& cfg,
                                 g2o::SparseOptimizer& open_loop_problem,
                                 const std::vector<g2o::EdgeSE2*>& loops);

class IPC
{
public :

    IPC(g2o::SparseOptimizer& open_loop_problem, const Config& cfg);
    ~IPC();

    bool agreementCheck(g2o::EdgeSE2* loop_candidate, struct Info& stats);

    inline const std::vector<std::pair<int, int>>& getMaxConsensusSet() const { return _max_consensus_set; }
    inline const std::vector<std::pair<int, int>>& getClusters() const { return _clusters; }

private :


    g2o::SparseOptimizer* _problem;

    std::vector<std::pair<int, int>> _max_consensus_set;
    std::vector<std::pair<int, int>> _clusters;
    std::vector<g2o::EdgeSE2*> _odom_edges;
    std::vector<g2o::OptimizableGraph::EdgeSet> _edgesxcl;
    std::vector<g2o::OptimizableGraph::EdgeSet> _edgesxcl_only_loops;

    double _fast_reject_th;
    int _fast_reject_iter_base;
    double _slow_reject_th;
    int _slow_reject_iter_base; 


};