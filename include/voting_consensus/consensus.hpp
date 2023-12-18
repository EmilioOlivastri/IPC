#pragma once

#include "voting_consensus/consensus_utils.hpp"
#include "voting_consensus/utils.hpp"

// Using full deterministic search 
void ipc(const Config& cfg,
         g2o::SparseOptimizer& open_loop_problem,
         const std::vector<g2o::EdgeSE2*>& loops,
         std::vector<g2o::EdgeSE2*>& inliers,
         std::vector<g2o::EdgeSE2*>& outliers,
         const std::string& output_file);