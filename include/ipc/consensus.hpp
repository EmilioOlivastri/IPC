#pragma once

#include "ipc/consensus_utils.hpp"

template <class EDGE, class VERTEX> class IPC
{
public :

    IPC(g2o::SparseOptimizer& open_loop_problem, const Config& cfg);
    ~IPC();

    bool agreementCheck(EDGE* loop_candidate);
    bool removeEdgeFromCnS(EDGE* edge);
    void addEdgeToCnS(EDGE* edge);

    const std::vector<EDGE*>& getMaxConsensusSet() const { return _max_consensus_set; }

private :

    std::pair<int, int> computeIndependentSubgraph(const EDGE& edge_candidate,
                                                   g2o::OptimizableGraph::EdgeSet& eset_independent);

    std::pair<int, int> getKBestLoops(const EDGE& edge_candidate, 
                                      g2o::OptimizableGraph::EdgeSet& kbest_edges);

    g2o::SparseOptimizer* _problem;

    std::vector<EDGE*> _max_consensus_set; 
    std::vector<EDGE*> _odom_edges;

    double _fast_reject_th;
    int _fast_reject_iter_base;
    double _slow_reject_th;
    int _slow_reject_iter_base;
    double _s_factor;

    bool _use_best_k_buddies;
    int _k_buddies;
};