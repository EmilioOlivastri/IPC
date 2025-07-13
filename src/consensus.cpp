#include "ipc/consensus.hpp"

using namespace std;
using namespace g2o;

// Main body of the outlier detection algorithm
template <class T, class EDGE, class VERTEX>
void simulating_incremental_data(const Config& cfg,
                                 SparseOptimizer& open_loop_problem,
                                 const vector<EDGE*>& loops)
{

    vector<T> poses_gt;
    readSolutionFile(poses_gt, cfg.ground_truth);
    double s_factor = cfg.s_factor;
    
    // 0 = outliers | 1 = inliers    
    int tot_hypothesis = loops.size();    
    vector<int> bucket(tot_hypothesis, 1);
    vector<pair<bool, EDGE*>> gt_loops;
    vector<EDGE*> sorted_loops;

    for (size_t idx = 0 ; idx < cfg.canonic_inliers; gt_loops.push_back(make_pair(true, loops[idx++])));
    for (size_t idx = cfg.canonic_inliers ; idx < loops.size(); gt_loops.push_back(make_pair(false, loops[idx++])));
    sort(gt_loops.begin(), gt_loops.end(), cmpTime);

    IPC<EDGE, VERTEX> ipc(open_loop_problem, cfg);

    // Simulating incremental addition of edges
    double avg_time = 0.0;
    cout << "Starting simulation of incremental dataset -> Displaying relative status : " << endl;
    cout << "S = " << s_factor << " | TH = " << cfg.fast_reject_th << endl;
    for ( size_t candidate_id = 0 ; candidate_id < tot_hypothesis ; ++candidate_id )
    {
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        bool consistent = ipc.agreementCheck(gt_loops[candidate_id].second);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        bucket[candidate_id] = consistent ? 1 : 0;
        sorted_loops.push_back(gt_loops[candidate_id].second);

        chrono::microseconds delta_time = chrono::duration_cast<chrono::microseconds>(end - begin);
        avg_time += delta_time.count() / 1000000.0;

        printProgress((double)(candidate_id + 1) / (double)tot_hypothesis);
    }
    cout << "\nCompleted!" << endl;

    vector<EDGE*> odom_edges;
    getProblemOdom<EDGE>(open_loop_problem, odom_edges);
    propagateGuess<EDGE, VERTEX>(open_loop_problem, 0, odom_edges.size(), odom_edges);

    OptimizableGraph::EdgeSet eset_gl;
    for ( size_t i = 0 ; i < odom_edges.size(); eset_gl.insert(odom_edges[i++]) )
        odom_edges[i]->setInformation(odom_edges[i]->information() / s_factor);
    for ( size_t i = 0 ; i < gt_loops.size(); ++i )
        if ( bucket[i] == 1 )
            eset_gl.insert(gt_loops[i].second);


    open_loop_problem.initializeOptimization(eset_gl);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    open_loop_problem.optimize(1000);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    chrono::microseconds delta_time = chrono::duration_cast<chrono::microseconds>(end - begin);
    double dt = delta_time.count() / 1000000.0;

    int tp = 0;
    int fp = 0;
    int tn = 0;
    int fn = 0;
    for ( size_t counter = 0 ; counter < gt_loops.size() ; ++counter )
        if ( gt_loops[counter].first && bucket[counter] == 1 ) ++tp;
        else if ( gt_loops[counter].first && bucket[counter] == 0 ) ++fn;
        else if ( !gt_loops[counter].first && bucket[counter] == 1 ) ++fp;
        else if ( !gt_loops[counter].first && bucket[counter] == 0 ) ++tn;

    float precision = tp / (float)(tp + fp);
    float recall = tp / (float)(tp + fn);

    const vector<pair<int, int>> clusters = ipc.getClusters();
    const vector<pair<int, int>> max_consensus_set = ipc.getMaxConsensusSet();

    for ( size_t counter = 0 ; counter < clusters.size() ; ++counter )    
        cout << "L("<< counter << ") = [" << clusters[counter].first << ", " << clusters[counter].second << "]\n";
    int size_max_consensus = 0;
    for ( size_t counter = 0 ; counter < gt_loops.size() ; size_max_consensus += bucket[counter++] );
    cout << "Size of MAX consistent set = " << size_max_consensus << endl;

    cout << "Avg Time x test = " << avg_time / tot_hypothesis << " [s]\n";
    cout << "Precision = " << precision << endl;        
    cout << "Recall = " << recall << endl;

    ofstream outfile;
    outfile.open(cfg.output.c_str()); 
    for (size_t it = 0; it < open_loop_problem.vertices().size(); ++it )
    {
        VERTEX* v = dynamic_cast<VERTEX*>(open_loop_problem.vertex(it));
        writeVertex(outfile, v);
    }
    outfile.close();


    string out2 = cfg.output.substr(0, cfg.output.size() - 3) + "PR";
    outfile.open(out2.c_str());
    outfile << precision << " " << recall << endl;
    outfile << avg_time << " " << avg_time / tot_hypothesis << endl;
    outfile.close();

    return;
}

template void simulating_incremental_data<Eigen::Isometry2d, EdgeSE2, VertexSE2>(const Config& cfg,
                                                                                 SparseOptimizer& open_loop_problem,
                                                                                 const vector<EdgeSE2*>& loops);
template void simulating_incremental_data<Eigen::Isometry3d, EdgeSE3, VertexSE3>(const Config& cfg,
                                                                                 SparseOptimizer& open_loop_problem,
                                                                                 const vector<EdgeSE3*>& loops);

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


template bool IPC<EdgeSE2, VertexSE2>::removeEdgeFromCnS(EdgeSE2* edge);
template bool IPC<EdgeSE3, VertexSE3>::removeEdgeFromCnS(EdgeSE3* edge);


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

template bool IPC<EdgeSE2, VertexSE2>::addEdgeToCnS(EdgeSE2* edge);
template bool IPC<EdgeSE3, VertexSE3>::addEdgeToCnS(EdgeSE3* edge);
