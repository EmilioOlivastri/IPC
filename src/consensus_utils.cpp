#include "ipc/consensus_utils.hpp"

using namespace std;
using namespace g2o;

template <class EDGE>
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
        if ( dynamic_cast<EDGE*>(e)->chi2() > th )
            return false;

    return true;
}

template bool isAgreeingWithCurrentState<EdgeSE2>(SparseOptimizer& problem, OptimizableGraph::EdgeSet& eset, double th, int iter_base);
template bool isAgreeingWithCurrentState<EdgeSE3>(SparseOptimizer& problem, OptimizableGraph::EdgeSet& eset, double th, int iter_base);

/*----------------------------------------------------------------*/

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

/*----------------------------------------------------------------*/

void fixAndFreeInternal(g2o::SparseOptimizer& problem, const std::set<int>& fixv, const std::set<int>& freev)
{
    for (auto& id : fixv)
        problem.vertex(id)->setFixed(true);

    for (auto& id : freev)
        problem.vertex(id)->setFixed(false);

    return;
}

/*----------------------------------------------------------------*/

template <class EDGE, class VERTEX>
void propagateCurrentGuess(SparseOptimizer& curr_estimate, int id_start, const vector<EDGE*>& odom)
{ 
    for ( size_t i = id_start + 1; i <= odom.size() ; ++i )
    {
        auto v1 = dynamic_cast<VERTEX*>(curr_estimate.vertex(i - 1));
        auto v2 = dynamic_cast<VERTEX*>(curr_estimate.vertex(i));
        v2->setEstimate(v1->estimate() * odom[i - 1]->measurement());
    }

    return;
}

template void propagateCurrentGuess<EdgeSE2, VertexSE2>(SparseOptimizer& curr_estimate, int id_start, const vector<EdgeSE2*>& odom);
template void propagateCurrentGuess<EdgeSE3, VertexSE3>(SparseOptimizer& curr_estimate, int id_start, const vector<EdgeSE3*>& odom);

/*----------------------------------------------------------------*/

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

/*----------------------------------------------------------------*/

template <class EDGE, class VERTEX>
void propagateGuess(SparseOptimizer& voting, int id1, int id2, const vector<EDGE*>& odom)
{ 
    // id1 is the gauge and origin of the voting system
    auto gauge = dynamic_cast<VERTEX*>(voting.vertex(id1));
    gauge->setFixed(true);
    gauge->setToOrigin();

    // Propagate guess for the vertices in the loop
    for ( size_t i = id1 + 1; i <= id2 ; ++i )
    {
        auto v1 = dynamic_cast<VERTEX*>(voting.vertex(i - 1));
        auto v2 = dynamic_cast<VERTEX*>(voting.vertex(i));
        v2->setFixed(false);
        v2->setEstimate(v1->estimate() * odom[i - 1]->measurement());
    }

    return;
}

template void propagateGuess<EdgeSE2, VertexSE2>(SparseOptimizer& voting, int id1, int id2, const vector<EdgeSE2*>& odom);
template void propagateGuess<EdgeSE3, VertexSE3>(SparseOptimizer& voting, int id1, int id2, const vector<EdgeSE3*>& odom);

/*----------------------------------------------------------------*/

template <class EDGE>
void robustifyVoters(int id1, int id2, double s_factor, vector<EDGE*>& voters)
{
    // Iterating through voters adding robust kernel
    for ( size_t j = id1 ; j < id2; ++j )
        voters[j]->setInformation(voters[j]->information() * s_factor);
    return;
}

template void robustifyVoters<EdgeSE2>(int id1, int id2, double s_factor, vector<EdgeSE2*>& voters);
template void robustifyVoters<EdgeSE3>(int id1, int id2, double s_factor, vector<EdgeSE3*>& voters);

/*----------------------------------------------------------------*/
template <class EDGE>
pair<int, int> getOrderedPair(const EDGE* edge)
{
    int v0_id = edge->vertices()[0]->id();
    int v1_id = edge->vertices()[1]->id();
    return make_pair(min(v0_id, v1_id), max(v0_id, v1_id));
}

template pair<int, int> getOrderedPair<EdgeSE2>(const EdgeSE2* edge);
template pair<int, int> getOrderedPair<EdgeSE3>(const EdgeSE3* edge);

/*----------------------------------------------------------------*/
double intersection_calc(const pair<int, int>& e1, const pair<int, int>& e2)
{
    // Start & End Idx of lp1
    int id1_st = e1.first;
    int id1_end = e1.second;

    // Start & End Idx of lp2
    int id2_st = e2.first;
    int id2_end = e2.second;

    // Compute intersection
    int intersection_st = id1_st > id2_st ? id1_st : id2_st;
    int intersection_end = id1_end < id2_end ? id1_end : id2_end;

    return static_cast<double>(intersection_end - intersection_st);
}

/*----------------------------------------------------------------*/
double union_calc(const pair<int, int>& e1, const pair<int, int>& e2)
{
    // Start & End Idx of lp1
    int id1_st = e1.first;
    int id1_end = e1.second;

    // Start & End Idx of lp2
    int id2_st = e2.first;
    int id2_end = e2.second;

    // Compute union
    int union_st = id1_st < id2_st ? id1_st : id2_st;
    int union_end = id1_end > id2_end ? id1_end : id2_end;

    return static_cast<double>(union_end - union_st);
}


/*----------------------------------------------------------------*/
double iou_calc(const pair<int, int>& e1, const pair<int, int>& e2)
{
    // Compute union
    double delta_union = union_calc(e1, e2);

    // Compute intersection
    double delta_intersection = intersection_calc(e1, e2);
    
    // Compute IoU (if negative no intersection)
    double iou = delta_intersection / delta_union;
    iou = iou < 0.0 ? 0.0 : iou;

    return iou;
}