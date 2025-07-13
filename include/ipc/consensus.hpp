#pragma once

#include "ipc/consensus_utils.hpp"

// Using full deterministic search
template <class T, class EDGE, class VERTEX>
void simulating_incremental_data(const Config& cfg,
                                 g2o::SparseOptimizer& open_loop_problem,
                                 const std::vector<EDGE*>& loops);

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> StateVector;

template <class EDGE, class VERTEX> class IPC
{
public :

    IPC(g2o::SparseOptimizer& open_loop_problem, const Config& cfg);
    ~IPC();

    bool agreementCheck(EDGE* loop_candidate);
    bool removeEdgeFromCnS(EDGE* edge);
    bool addEdgeToCnS(EDGE* edge);

    const std::vector<std::pair<int, int>>& getMaxConsensusSet() const { return _max_consensus_set; }
    const std::vector<std::pair<int, int>>& getClusters() const { return _clusters; }

private :

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
    double _s_factor;
};