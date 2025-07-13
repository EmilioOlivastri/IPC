#pragma once

#include "ipc/consensus.hpp"

template <class T, class EDGE, class VERTEX>
void simulating_incremental_data(const Config& cfg,
                                 g2o::SparseOptimizer& open_loop_problem,
                                 const std::vector<EDGE*>& loops);

                                 