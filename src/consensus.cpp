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

    _use_best_k_buddies = cfg.use_best_k_buddies;
    _k_buddies = _use_best_k_buddies ? cfg.k_buddies : -1;
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
    pair<int, int> clust_ext =  _use_best_k_buddies ? getKBestLoops(*loop_candidate, eset_independent) :
                                computeIndependentSubgraph(*loop_candidate, eset_independent);

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
    pair<int, int> cand = getOrderedPair(&edge_candidate);

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

            const auto& trusted_edge = getOrderedPair(_max_consensus_set[edge_cidx]);

            // Compute intersection with the current extremes
            int intersection = static_cast<int>(intersection_calc(extremes, trusted_edge));

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


template <class EDGE, class VERTEX> 
pair<int, int> IPC<EDGE, VERTEX>::getKBestLoops(const EDGE& edge_candidate, 
                                                OptimizableGraph::EdgeSet& kbest_edges)
{
    // Assuring the output set is empty
    kbest_edges.clear();

    // Get the vertices ids of the candidate edge
    pair<int, int> cand = getOrderedPair(&edge_candidate);

    // Initialize the extremes of the independent subgraph using the candidate edge
    pair<int, int> extremes(cand.first, cand.second);

    // Compute IoU for all the voters
    int iters = min(_k_buddies, (int)_max_consensus_set.size());
    vector<bool> included_edges(_max_consensus_set.size(), false);
    for (int idx = 0; idx < iters; ++idx)
    {
        double max_iou = -1.0;
        double min_union = std::numeric_limits<double>::max();
        int max_iou_idx = -1;
        for (size_t vt_loop_idx = 0; vt_loop_idx < _max_consensus_set.size(); ++vt_loop_idx)
        {
            // Skip if already included
            if (included_edges[vt_loop_idx]) continue;

            pair<int, int> vt_pair = getOrderedPair(_max_consensus_set[vt_loop_idx]);

            // Compute IoU
            double iou_score = iou_calc(extremes, vt_pair);

            // Skip if IoU is STRICTLY LESS than the maximum found so far
            if (iou_score < max_iou || iou_score == 0.0) continue;

            // Compute the union of the vertices ids, we prefer the one with the smallest union
            double union_score = union_calc(extremes, vt_pair);

            // Skip if the iou is the same but the union is not better than the current best
            /**/
            if (iou_score ==  max_iou && union_score >= min_union) 
                continue;  
            /**/

            // Update maximum IoU and index
            max_iou = iou_score;
            max_iou_idx = vt_loop_idx;
            min_union = union_score;
        }

        if (max_iou_idx < 0 ) continue;

        included_edges[max_iou_idx] = true;
        kbest_edges.insert(_max_consensus_set[max_iou_idx]);

        // Update extremes
        pair<int, int> vt_pair = getOrderedPair(_max_consensus_set[max_iou_idx]);
        extremes.first = min(extremes.first, vt_pair.first);
        extremes.second = max(extremes.second, vt_pair.second);
    }

    cout << "K-Best Loops: " << kbest_edges.size() << endl;
    cout << "For the candidate edge: [" << cand.first << ", " << cand.second << "]\n";
    cout << "This are the best loops:\n";
    for (const auto& hyper_edge : kbest_edges)
    {
        auto edge = dynamic_cast<const EDGE*>(hyper_edge);
        pair<int, int> loop_pair = getOrderedPair(edge);
        cout << "L(" << loop_pair.first << ", " << loop_pair.second << ")\n";
    }
    std::cout << "Extremes: [" << extremes.first << ", " << extremes.second << "]\n";
    cout << "----------------------------------------\n";

    return extremes;    
}


template class IPC<EdgeSE2, VertexSE2>;
template class IPC<EdgeSE3, VertexSE3>;