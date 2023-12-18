#include "ipc/consensus_utils.hpp"

using namespace std;
using namespace g2o;

bool isAgreeingWithCurrentState(SparseOptimizer& problem, OptimizableGraph::EdgeSet& eset, double th, int iter_base)
{
    // Generate map hypothesis closed-loop subproblem
    problem.initializeOptimization(eset);
    problem.computeActiveErrors();
    int iter = iter_base;
    iter = eset.size() > 100 ? iter * 5 : iter ;
    problem.optimize(iter);
    problem.computeActiveErrors();

    for (auto& e : eset )
        if ( dynamic_cast<EdgeSE2*>(e)->chi2() > th )
            return false;

    return true;
}

void fixComplementary(SparseOptimizer& problem, int start_id, int end_id)
{
    // Fixing vertices before
    for (size_t id = 0; id <= start_id ; ++id )
        problem.vertex(id)->setFixed(true);

    for (size_t id = start_id + 1; id <= end_id ; ++id )
        problem.vertex(id)->setFixed(false);

    // Fixing all the vertices after
    for (size_t id = end_id + 1; id < problem.vertices().size() ; ++id )
        problem.vertex(id)->setFixed(true);

    return;
}

void propagateCurrentGuess(SparseOptimizer& curr_estimate, int id_start, const vector<EdgeSE2*>& odom)
{ 
    for ( size_t i = id_start + 1; i <= odom.size() ; ++i )
    {
        auto v1 = dynamic_cast<VertexSE2*>(curr_estimate.vertex(i - 1));
        auto v2 = dynamic_cast<VertexSE2*>(curr_estimate.vertex(i));
        v2->setEstimate(v1->estimate() * odom[i - 1]->measurement());

    }

    return;
}

vector<int> getClusterId(const vector<pair<int, int>>& clusters, const pair<int, int>& candidate)
{
    int intersection = -1;
    vector<int> clusters_id;
    for ( size_t cl_id = 0 ; cl_id < clusters.size(); ++cl_id )
    {
        // Check if out of the extremes of the cluster
        bool out_right = clusters[cl_id].second <= candidate.first ? true : false;
        bool out_left = clusters[cl_id].first >= candidate.second ? true : false;

        if ( out_right || out_left ) continue;

        clusters_id.push_back(cl_id);
    }

    return clusters_id;
}

void propagateGuess(SparseOptimizer& voting, int id1, int id2, const vector<EdgeSE2*>& odom)
{ 
    // id1 is the gauge and origin of the voting system
    auto gauge = dynamic_cast<VertexSE2*>(voting.vertex(id1));
    gauge->setFixed(true);
    gauge->setEstimate(SE2(0, 0, 0));

    // Propagate guess for the vertices in the loop
    for ( size_t i = id1 + 1; i <= id2 ; ++i )
    {
        auto v1 = dynamic_cast<VertexSE2*>(voting.vertex(i - 1));
        auto v2 = dynamic_cast<VertexSE2*>(voting.vertex(i));
        v2->setFixed(false);
        v2->setEstimate(v1->estimate() * odom[i - 1]->measurement());
    }

    return;
}

void robustifyVoters(int id1, int id2, double sqrt_th, vector<EdgeSE2*>& voters)
{
    // Iterating through voters adding robust kernel
    for ( size_t j = id1 ; j < id2; ++j )
        voters[j]->setInformation(voters[j]->information() * 10);
    return;
}