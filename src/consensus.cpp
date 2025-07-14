#include "ipc/consensus.hpp"

using namespace std;
using namespace g2o;


// Constructor
template <class EDGE, class VERTEX> 
IPC<EDGE, VERTEX>::IPC(g2o::SparseOptimizer& open_loop_problem, const Config& cfg)
{
    // Inizialization
    _problem = &open_loop_problem;
    getProblemOdom<EDGE>(*_problem, _odom_edges);

    std::sort(_odom_edges.begin(), _odom_edges.end(), cmpEdgesID);

    // Robustify simple voters --> only once
    _s_factor = cfg.s_factor;
    double th = cfg.fast_reject_th;
    double sqrt_th = sqrt(th);
    robustifyVoters<EDGE>(0, _odom_edges.size(), _s_factor, _odom_edges);

    propagateGuess<EDGE, VERTEX>(*_problem, 0, _odom_edges.size(), _odom_edges);
    
    store<VERTEX>(*_problem);

    _max_consensus_set.clear();

    _fast_reject_th = cfg.fast_reject_th;
    _fast_reject_iter_base = cfg.fast_reject_iter_base;
    _slow_reject_th = cfg.slow_reject_th;
    _slow_reject_iter_base = cfg.slow_reject_iter_base;
}

template <class EDGE, class VERTEX> 
IPC<EDGE, VERTEX>::~IPC()
{
    _problem->clear();
    _max_consensus_set.clear();
}

template <class EDGE, class VERTEX>
bool IPC<EDGE, VERTEX>::agreementCheck(EDGE* loop_candidate)
{
    // Estimate the independent subgraph given the loop candidate
    OptimizableGraph::EdgeSet eset_independent;
    pair<int, int> clust_ext = computeIndependentSubgraph(*loop_candidate, eset_independent);

    // eset_independent.size() == 0 --> no intersection found
    bool intersection_found = (eset_independent.size() > 0);
    double reject_th = intersection_found ? _slow_reject_th : _fast_reject_th;
    int reject_iter_base = intersection_found ? _slow_reject_iter_base : _fast_reject_iter_base;
    
    // Adding the reliable voters
    for ( size_t j = clust_ext.first ; j < clust_ext.second; eset_independent.insert(_odom_edges[j++]) );
    eset_independent.insert(loop_candidate);

    // Generate map hypothesis closed-loop subproblem
    store<VERTEX>(*_problem);
    fixComplementary(*_problem, clust_ext.first, clust_ext.second);

    bool is_agreeing = isAgreeingWithCurrentState<EDGE>(*_problem, eset_independent, reject_th, reject_iter_base);
    if (!is_agreeing)
    {
        restore<VERTEX>(*_problem);
        return false;
    }

    discard<VERTEX>(*_problem);
    _max_consensus_set.push_back(loop_candidate);
    propagateCurrentGuess<EDGE, VERTEX>(*_problem, clust_ext.second, _odom_edges);

    return true;
        
}

template <class EDGE, class VERTEX>
bool IPC<EDGE, VERTEX>::removeEdgeFromCnS(EDGE* edge)
{
    // Check if the edge is in the consensus set
    int v0_id = min(edge->vertices()[0]->id(), edge->vertices()[1]->id());
    int v1_id = max(edge->vertices()[0]->id(), edge->vertices()[1]->id());
    for (auto it = _max_consensus_set.begin(); it != _max_consensus_set.end(); ++it)
    {
        EDGE* trusty_edge = *it;
        int vt0_id = min(trusty_edge->vertices()[0]->id(), trusty_edge->vertices()[1]->id());
        int vt1_id = max(trusty_edge->vertices()[0]->id(), trusty_edge->vertices()[1]->id());

        if (v0_id != vt0_id || v1_id != vt1_id )
            continue;

        _max_consensus_set.erase(it);
        return true;
    }

    return false;
}

template <class EDGE, class VERTEX>
void IPC<EDGE, VERTEX>::addEdgeToCnS(EDGE* edge)
{
    // Check if the edge is already in the consensus set
    int v0_id = min(edge->vertices()[0]->id(), edge->vertices()[1]->id());
    int v1_id = max(edge->vertices()[0]->id(), edge->vertices()[1]->id());
    for (const auto& trusty_edge : _max_consensus_set)
    {
        int vt0_id = min(trusty_edge->vertices()[0]->id(), 
                         trusty_edge->vertices()[1]->id());
        int vt1_id = max(trusty_edge->vertices()[0]->id(), 
                         trusty_edge->vertices()[1]->id());

        if (v0_id == vt0_id && v1_id == vt1_id) return;
    }
    
    _max_consensus_set.push_back(edge);

    // Sort the consensus set by vertex ids
    std::sort(_max_consensus_set.begin(), _max_consensus_set.end(), cmpEdgesTime);

    return;
}

template <class EDGE, class VERTEX>
pair<int, int> IPC<EDGE, VERTEX>::computeIndependentSubgraph(const EDGE& edge_candidate,
                                                             OptimizableGraph::EdgeSet& eset_independent)
{
    // Get the vertices ids of the candidate edge
    int v0_id = min(edge_candidate.vertices()[0]->id(), 
                    edge_candidate.vertices()[1]->id());
    int v1_id = max(edge_candidate.vertices()[0]->id(), 
                    edge_candidate.vertices()[1]->id());
    pair<int, int> cand(v0_id, v1_id);

    // Initialize the extremes of the independent subgraph using the candidate edge
    pair<int, int> extremes(cand.first, cand.second);

    // Iterate through the consensus set to check for edge intersections
    bool found_new_extremes = true;
    vector<bool> included_edges(_max_consensus_set.size(), false);
    while (found_new_extremes)
    {
        // Check for intersections with the current extremes and skip already included edges
        found_new_extremes = false;

        for (size_t edge_cidx = 0; edge_cidx < _max_consensus_set.size(); ++edge_cidx)
        {
            // If the edge is already included, skip it
            if (included_edges[edge_cidx]) continue;

            int vt0_id = min(_max_consensus_set[edge_cidx]->vertices()[0]->id(),
                             _max_consensus_set[edge_cidx]->vertices()[1]->id());
            int vt1_id = max(_max_consensus_set[edge_cidx]->vertices()[0]->id(),
                             _max_consensus_set[edge_cidx]->vertices()[1]->id());
            const auto& trusted_edge = make_pair(vt0_id, vt1_id);

            // Compute intersection with the current extremes
            int intersection = min(trusted_edge.second, extremes.second) - max(trusted_edge.first, extremes.first);
            
            if (intersection <= 0) continue;  // No intersection

            // Update extremes if an intersection is found
            extremes.first = min(extremes.first, trusted_edge.first);
            extremes.second = max(extremes.second, trusted_edge.second);
            included_edges[edge_cidx] = true;
            found_new_extremes = true;
            eset_independent.insert(_max_consensus_set[edge_cidx]);
        }
    }

    return extremes;
}


template class IPC<EdgeSE2, VertexSE2>;
template class IPC<EdgeSE3, VertexSE3>;