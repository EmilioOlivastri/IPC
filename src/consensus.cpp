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
    _clusters.clear();
    _edgesxcl.clear();
    _edgesxcl_only_loops.clear();

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
    _clusters.clear();
    _edgesxcl.clear();
    _edgesxcl_only_loops.clear();
}

template <class EDGE, class VERTEX>
bool IPC<EDGE, VERTEX>::agreementCheck(EDGE* loop_candidate)
{
    // Start & End Idx of Simple Loop
    int v1_id = loop_candidate->vertices()[0]->id();
    int v2_id = loop_candidate->vertices()[1]->id();
    int id1 = min(v1_id, v2_id);
    int id2 = max(v1_id, v2_id);

    pair<int, int> candidate(id1, id2);
    vector<int> vcl_id = getClusterId(_clusters, candidate);
    
    // No intersection found, creation of new cluster
    if ( vcl_id.size() == 0 )
    {
        // Adding the reliable voters
        OptimizableGraph::EdgeSet eset, eset_loop;
        for ( size_t j = id1 ; j < id2; eset.insert(_odom_edges[j++]) );
        eset.insert(loop_candidate);

        // Generate map hypothesis closed-loop subproblem
        store<VERTEX>(*_problem);
        fixComplementary(*_problem, id1, id2);

        if (!isAgreeingWithCurrentState<EDGE>(*_problem, eset, _fast_reject_th, _fast_reject_iter_base))
        {
            restore<VERTEX>(*_problem);
            return false;
        }
        else
        {
            discard<VERTEX>(*_problem);
            _clusters.push_back(candidate);
            _edgesxcl.push_back(eset);
            eset_loop.insert(loop_candidate);
            _edgesxcl_only_loops.push_back(eset_loop);
            _max_consensus_set.push_back(candidate);
            propagateCurrentGuess<EDGE, VERTEX>(*_problem, id2, _odom_edges);

            return true;
        }
    }
        
    // If the loop connects more clusters, then use the measurments from those clusters
    OptimizableGraph::EdgeSet eset, eset_loops, test;

    // Create a fake vertex for sorting the clusters & fusing cluster edges
    vector<pair<int, int>> sorted_intersect_clust;
    for ( size_t idx = 0; idx < vcl_id.size(); ++idx )
    {
        sorted_intersect_clust.push_back(_clusters[vcl_id[idx]]);
        eset.insert(_edgesxcl[vcl_id[idx]].begin(), _edgesxcl[vcl_id[idx]].end());
        eset_loops.insert(_edgesxcl_only_loops[vcl_id[idx]].begin(), _edgesxcl_only_loops[vcl_id[idx]].end());
    }
    sort(sorted_intersect_clust.begin(), sorted_intersect_clust.end(), cmpSecond);

    // Adding odom edges between connected clusters
    for ( size_t idx = 0; idx < sorted_intersect_clust.size() - 1; ++idx )
        for ( size_t j = sorted_intersect_clust[idx].second ; j < sorted_intersect_clust[idx + 1].first; eset.insert(_odom_edges[j++]) );

    // New cluster beggining/ending
    int new_cl_start = sorted_intersect_clust[0].first;
    int new_cl_end = sorted_intersect_clust[sorted_intersect_clust.size() - 1].second;

    // Adding extra odom edges to the left and to the right
    for ( size_t idx = candidate.first ; idx < new_cl_start; eset.insert(_odom_edges[idx++]) );
    for ( size_t idx = new_cl_end ; idx < candidate.second ; eset.insert(_odom_edges[idx++]) );

    // Compute new cluster extremes
    new_cl_start = new_cl_start < candidate.first ? new_cl_start : candidate.first;
    new_cl_end = new_cl_end > candidate.second ? new_cl_end : candidate.second;

    // Finally! Adding new_candidate edge!!
    eset.insert(loop_candidate);
    if ( vcl_id.size() > 1 )
        propagateGuess<EDGE, VERTEX>(*_problem, 0, _odom_edges.size(), _odom_edges);

    // Generate map hypothesis closed-loop subproblem
    store<VERTEX>(*_problem);

    int start = new_cl_start;
    int end = new_cl_end;
    test = eset;

    fixComplementary(*_problem, start, end);
    eset_loops.insert(loop_candidate);

    if (!isAgreeingWithCurrentState<EDGE>(*_problem, test, 
                                           _slow_reject_th, 
                                           _slow_reject_iter_base))
    {
        restore<VERTEX>(*_problem);
        return false;
    }
    
    discard<VERTEX>(*_problem);            
    _max_consensus_set.push_back(candidate);
    propagateCurrentGuess<EDGE, VERTEX>(*_problem, start, _odom_edges);

    if ( vcl_id.size() > 1 )
    {
        vector<pair<int, int>> clusters_tmp;
        vector<OptimizableGraph::EdgeSet> edgesxcl_tmp;
        vector<OptimizableGraph::EdgeSet> edgesxcl_only_loops_tmp;

        for ( size_t cidx = 0; cidx < _clusters.size(); ++cidx )
        {
            // Check if out of the extremes of the cluster
            bool out_right = _clusters[cidx].second <= new_cl_start ? true : false;
            bool out_left = _clusters[cidx].first >= new_cl_end ? true : false;

            if ( !(out_right || out_left) ) continue;

            // Adding if independent
            clusters_tmp.push_back(_clusters[cidx]);
            edgesxcl_tmp.push_back(_edgesxcl[cidx]);
            edgesxcl_only_loops_tmp.push_back(_edgesxcl_only_loops[cidx]);
        }
        _clusters = clusters_tmp;
        _edgesxcl = edgesxcl_tmp;
        _edgesxcl_only_loops = edgesxcl_only_loops_tmp;

        _clusters.push_back(make_pair(new_cl_start, new_cl_end));
        _edgesxcl.push_back(eset);
        _edgesxcl_only_loops.push_back(eset_loops);
    }
    else 
    {
        _clusters[vcl_id[0]].first = new_cl_start;
        _clusters[vcl_id[0]].second = new_cl_end;
        _edgesxcl[vcl_id[0]] = eset;
        _edgesxcl_only_loops[vcl_id[0]] = eset_loops;
    }

    return true;
        
}

template <class EDGE, class VERTEX>
bool IPC<EDGE, VERTEX>::removeEdgeFromCnS(EDGE* edge)
{
    // FIND CLUSTER OF BELONGING
    int v1_id = edge->vertices()[0]->id();
    int v2_id = edge->vertices()[1]->id();
    int id_src = min(v1_id, v2_id);
    int id_dst = max(v1_id, v2_id);

    pair<int, int> candidate(id_src, id_dst);
    vector<int> vcl_id = getClusterId(_clusters, candidate);

    // When removing an edge, the cluster is only 1 for sure
    int clust_id = vcl_id[0];
    pair<int, int> clust_ext = _clusters[clust_id];
    _edgesxcl_only_loops[clust_id].erase(edge);
    _edgesxcl[clust_id].erase(edge);

    /**/
    if ( _edgesxcl[clust_id].size() == 0 )
    {
        _clusters.erase(_clusters.begin() + clust_id);
        _edgesxcl.erase(_edgesxcl.begin() + clust_id);
        _edgesxcl_only_loops.erase(_edgesxcl_only_loops.begin() + clust_id);
        return true;
    }

    // Estimate new cluster extremes
    int new_cl_start = INT_MAX;
    int new_cl_end = -1;
    for (const auto& el : _edgesxcl_only_loops[clust_id])
    {
        int v1_id = el->vertices()[0]->id();
        int v2_id = el->vertices()[1]->id();
        new_cl_start = min(new_cl_start, min(v1_id, v2_id));
        new_cl_end = max(new_cl_end, max(v1_id, v2_id));
    }
    _clusters[clust_id].first = new_cl_start;
    _clusters[clust_id].second = new_cl_end;

    // Removing extra odom edges to the left and to the right
    if ( new_cl_start > clust_ext.first )
        for ( size_t idx = clust_ext.first ; idx < new_cl_start; _edgesxcl[clust_id].erase(_odom_edges[idx++]) );

    if ( new_cl_end < clust_ext.second )
        for ( size_t idx = new_cl_end ; idx < clust_ext.second; _edgesxcl[clust_id].erase(_odom_edges[idx++]) );

    return true;
}

template <class EDGE, class VERTEX>
bool IPC<EDGE, VERTEX>::addEdgeToCnS(EDGE* edge)
{
    // FIND CLUSTER OF BELONGING
    int v1_id = edge->vertices()[0]->id();
    int v2_id = edge->vertices()[1]->id();
    int id_src = min(v1_id, v2_id);
    int id_dst = max(v1_id, v2_id);

    pair<int, int> candidate(id_src, id_dst);
    vector<int> vcl_id = getClusterId(_clusters, candidate);

    OptimizableGraph::EdgeSet eset, eset_loops;

    // Case it needs a new cluster (unlikely)
    if (vcl_id.size() == 0)
    {
        for ( size_t j = id_src ; j < id_dst; eset.insert(_odom_edges[j++]) );
        eset.insert(edge);

        _clusters.push_back(candidate);
        _edgesxcl.push_back(eset);
        eset_loops.insert(edge);
        _edgesxcl_only_loops.push_back(eset_loops);
        
        return true;
    }

    if ( vcl_id.size() == 1 )
    {
        // Adding extra odom edges to the left and to the right
        int clust_id = vcl_id[0];
        pair<int, int> clust_ext = _clusters[clust_id];

        for ( size_t idx = candidate.first ; idx < clust_ext.first; _edgesxcl[clust_id].insert(_odom_edges[idx++]) );
        for ( size_t idx = clust_ext.second ; idx < candidate.second ; _edgesxcl[clust_id].insert(_odom_edges[idx++]) );

        _clusters[clust_id].first = candidate.first < clust_ext.first ? candidate.first : clust_ext.first;
        _clusters[clust_id].second = candidate.second > clust_ext.second ? candidate.second : clust_ext.second;
        _edgesxcl_only_loops[clust_id].insert(edge);
        _edgesxcl[clust_id].insert(edge);
        _max_consensus_set.push_back(candidate);

        return true;
    }

    /**/
    if ( vcl_id.size() > 1 )
    {
        // Case it needs to add to an existing cluster
        // If the loop connects more clusters, then use the measurments from those clusters
        vector<pair<int, int>> sorted_intersect_clust;
        for ( size_t idx = 0; idx < vcl_id.size(); ++idx )
        {
            sorted_intersect_clust.push_back(_clusters[vcl_id[idx]]);
            eset.insert(_edgesxcl[vcl_id[idx]].begin(), _edgesxcl[vcl_id[idx]].end());
            eset_loops.insert(_edgesxcl_only_loops[vcl_id[idx]].begin(), _edgesxcl_only_loops[vcl_id[idx]].end());
        }
        sort(sorted_intersect_clust.begin(), sorted_intersect_clust.end(), cmpSecond);

        // Adding odom edges between connected clusters
        for ( size_t idx = 0; idx < sorted_intersect_clust.size() - 1; ++idx )
            for ( size_t j = sorted_intersect_clust[idx].second ; j < sorted_intersect_clust[idx + 1].first; eset.insert(_odom_edges[j++]) );

        // New cluster beggining/ending
        int new_cl_start = sorted_intersect_clust[0].first;
        int new_cl_end = sorted_intersect_clust[sorted_intersect_clust.size() - 1].second;

        for ( size_t idx = candidate.first ; idx < new_cl_start; eset.insert(_odom_edges[idx++]) );
        for ( size_t idx = new_cl_end ; idx < candidate.second ; eset.insert(_odom_edges[idx++]) );

        // Compute new cluster extremes
        new_cl_start = new_cl_start < candidate.first ? new_cl_start : candidate.first;
        new_cl_end = new_cl_end > candidate.second ? new_cl_end : candidate.second;

        // Finally! Adding new_candidate edge!!
        eset.insert(edge);
        eset_loops.insert(edge);
        _max_consensus_set.push_back(candidate);

        vector<pair<int, int>> clusters_tmp;
        vector<OptimizableGraph::EdgeSet> edgesxcl_tmp;
        vector<OptimizableGraph::EdgeSet> edgesxcl_only_loops_tmp;

        for ( size_t cidx = 0; cidx < _clusters.size(); ++cidx )
        {
            // Check if out of the extremes of the cluster
            bool out_right = _clusters[cidx].second <= new_cl_start ? true : false;
            bool out_left = _clusters[cidx].first >= new_cl_end ? true : false;

            if ( !(out_right || out_left) ) continue;

            // Adding if independent
            clusters_tmp.push_back(_clusters[cidx]);
            edgesxcl_tmp.push_back(_edgesxcl[cidx]);
            edgesxcl_only_loops_tmp.push_back(_edgesxcl_only_loops[cidx]);
        }
        _clusters = clusters_tmp;
        _edgesxcl = edgesxcl_tmp;
        _edgesxcl_only_loops = edgesxcl_only_loops_tmp;

        _clusters.push_back(make_pair(new_cl_start, new_cl_end));
        _edgesxcl.push_back(eset);
        _edgesxcl_only_loops.push_back(eset_loops);
    }
    /**/

    return false;
}


template class IPC<EdgeSE2, VertexSE2>;
template class IPC<EdgeSE3, VertexSE3>;