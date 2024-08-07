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
        odom_edges[i]->setInformation(odom_edges[i]->information() * 0.1);
    for ( size_t i = 0 ; i < gt_loops.size(); ++i )
        if ( bucket[i] == 1 )
            eset_gl.insert(gt_loops[i].second);

    if ( cfg.use_recovery )
    {
        cout << "Starting Recovery Procedure" << endl;
        ipc.recoverySC(bucket, sorted_loops);
        cout << "Recovery Procedure Completed" << endl;
    }

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
    cout << "Avg Time x test = " << avg_time / tot_hypothesis << " [s]\n";
    cout << "Precision = " << precision << endl;        
    cout << "Recall = " << recall << endl;        
    cout << "Size of MAX consistent set = " << max_consensus_set.size() << endl;

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
    outfile << dt << endl;
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

    // Robustify simple voters --> only once
    double th = cfg.fast_reject_th;
    double sqrt_th = sqrt(th);
    robustifyVoters<EDGE>(0, _odom_edges.size(), sqrt_th, _odom_edges);

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
    _use_best_k_buddies = cfg.use_best_k_buddies;
    _k_buddies = _use_best_k_buddies ? cfg.k_buddies : -1;

    _use_recovery = true;
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
    int id1 = loop_candidate->vertices()[0]->id();
    int id2 = loop_candidate->vertices()[1]->id();
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
    if ( _use_best_k_buddies )
    {
        // Compute the top K buddies
        pair<int, int> clust_buddy = getTopKVoters(candidate, eset_loops, _k_buddies, test);
        
        // Just voting set edges, adding the odometry ones
        for ( size_t idx = clust_buddy.first; idx < clust_buddy.second; test.insert(_odom_edges[idx++]));

        // Compute new cluster extremes
        start = clust_buddy.first;
        end = clust_buddy.second;

        test.insert(loop_candidate);
    }
    else test = eset;

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
double IPC<EDGE, VERTEX>::iou(const pair<int, int>& lp1, const pair<int, int>& lp2)
{
    // Start & End Idx of lp1
    int id1_st = lp1.first;
    int id1_end = lp1.second;

    // Start & End Idx of lp2
    int id2_st = lp2.first;
    int id2_end =lp2.second;

    // Compute union
    int union_st = id1_st < id2_st ? id1_st : id2_st;
    int union_end = id1_end > id2_end ? id1_end : id2_end;

    // Compute intersection
    int intersection_st = id1_st > id2_st ? id1_st : id2_st;
    int intersection_end = id1_end < id2_end ? id1_end : id2_end;

    // Compute IoU (if negative no intersection)
    double iou = (double)(intersection_end - intersection_st) / (double)(union_end - union_st);
    iou = iou < 0.0 ? 0.0 : iou;

    return iou;
}


template <class EDGE, class VERTEX> 
pair<int, int> IPC<EDGE, VERTEX>::getTopKVoters(const pair<int, int>& candidate, 
                                                const OptimizableGraph::EdgeSet& eset_loops,
                                                int k, OptimizableGraph::EdgeSet& voters)
{
    // Assuring the the voters set is empty
    voters.clear();

    // Compute IoU for all the voters
    vector<pair<double, OptimizableGraph::Edge*>> score_vec;
    for ( const auto& vt_loop : eset_loops )
    {
        // Compute IoU
        pair<int, int> vt_pair(vt_loop->vertices()[0]->id(), 
                               vt_loop->vertices()[1]->id());

        double score = iou(candidate, vt_pair);
        score_vec.push_back(make_pair(score, dynamic_cast<OptimizableGraph::Edge*>(vt_loop)));
    }

    // Extract the top K voters or the first ones that have IoU > 0.0
    sort(score_vec.begin(), score_vec.end(), cmpScores);
    int min_k = min(k, (int)score_vec.size());
    
    pair<int, int> cluster = candidate;
    for ( int idx = 0 ; idx < min_k && score_vec[idx].first > 0.0; ++idx )
    {
        voters.insert(score_vec[idx].second);
        cluster.first = min(cluster.first, score_vec[idx].second->vertices()[0]->id());
        cluster.second = max(cluster.second, score_vec[idx].second->vertices()[1]->id());
    }

    return cluster;    
}