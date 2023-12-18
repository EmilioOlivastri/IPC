#pragma once

#include "voting_consensus/utils.hpp"

/*
Method that verifies if inliers accept the new measurment. It is based Chi-Squared Test.
@ param problem : sub_problem that needs to be optimized;
@ param eset : edges that partecipate to the optimized;
@ param th : threshold for the Chi-Squared Test;
@ return val : true if agree, false if not agree;
*/
bool isAgreeingWithCurrentState(g2o::SparseOptimizer& problem, g2o::OptimizableGraph::EdgeSet& eset, double th, int iter_base, int& eff_iters);

/*
Method that fixes the vertices that do not participate to the optimization.
@ param problem  : problem that needs to be solved/optimized;
@ param start_id : first vertex that partecipates to the opt;
@ param end_id   : last vertex that partecipates to the opt;
*/
void fixComplementary(g2o::SparseOptimizer& problem, int start_id, int end_id);


/*
Method that propagates the current estimate to the rest of the problem.
@ param curr_estimate  : current problem solution;
@ param id_start       : index of the last vertex of the active subgraph;
@ param odom           : odometry edges that will be used for propagation;
*/
void propagateCurrentGuess(g2o::SparseOptimizer& curr_estimate, int id_start, const std::vector<g2o::EdgeSE2*>& odom);

/*
Method that returns which cluster the candidate is intersecting.
@ param clusters : sub_problem that needs to be optimized;
@ param candidate : edges that partecipate to the optimized;
@ return : vector of clusters' id that intersect the candidate;
*/
std::vector<int> getClusterId(const std::vector<std::pair<int, int>>& clusters, const std::pair<int, int>& candidate);

/*
Propagation of the guess obtained from the odometry. 
It sets the vertex with id = id1 to the origin.
@ param graph : representation of the full problem;
@ param id1 : starting vertex id of the subgraph;
@ param id2 : ending vertex id of the subgraph;
@ param odom : odometry edges that will be used for propagation;
*/
void propagateGuess(g2o::SparseOptimizer& graph, int id1, int id2, const std::vector<g2o::EdgeSE2*>& odom);

/*
Method increases the weight of the odometry edges partecipating to the optimization.
@ param id1 : starting vertex id of the subgraph;
@ param id2 : ending vertex id of the subgraph;
@ param sqrt_th : threshold for the Chi-Squared Test;
@ param voters : edges that partecipate to the optimization;
@ return val : true if agree, false if not agree;
*/
void robustifyVoters(int id1, int id2, double sqrt_th, std::vector<g2o::EdgeSE2*>& voters);