#include "voting_consensus/consensus.hpp"

using namespace std;
using namespace g2o;

// Main body of the outlier detection algorithm
void ipc(const Config& cfg,
         SparseOptimizer& open_loop_problem,
         const vector<EdgeSE2*>& loops,
         vector<EdgeSE2*>& inliers,
         vector<EdgeSE2*>& outliers,
         const string& output_file)
{
    std::vector<Eigen::Isometry2d> poses_gt;
    readSolutionFile(poses_gt, cfg.ground_truth);
    
    // 0 = outliers | 1 = not tested | > 1 = size of consensus set
    int tot_hypothesis = loops.size();    
    vector<int> bucket(tot_hypothesis, 1);
    vector<pair<bool, EdgeSE2*>> gt_loops;

    for (size_t idx = 0 ; idx < cfg.canonic_inliers; gt_loops.push_back(make_pair(true, loops[idx++])));
    for (size_t idx = cfg.canonic_inliers ; idx < loops.size(); gt_loops.push_back(make_pair(false, loops[idx++])));
    sort(gt_loops.begin(), gt_loops.end(), cmpTime);

    // Extract the odom edges
    // Questa funzione extra si puo' evitare
    // se diamo in input direttamente gli odom edges
    vector<EdgeSE2*> odom_edges;
    getProblemOdom(open_loop_problem, odom_edges);

    // Robustify simple voters --> only once
    double th = cfg.fast_reject_th;
    double sqrt_th = sqrt(th);
    robustifyVoters(0, odom_edges.size(), sqrt_th, odom_edges);

    vector<pair<int, int>> max_consensus_set;
    vector<pair<int, int>> clusters;
    vector<OptimizableGraph::EdgeSet> edgesxcl;
    vector<OptimizableGraph::EdgeSet> edgesxcl_only_loops;
    propagateGuess(open_loop_problem, 0, odom_edges.size(), odom_edges);
    store(open_loop_problem);

    // STUDY DISTRIBUTION OF REJECTIONS / OUTLIERS IN TIME
    vector<struct Info> statistics(tot_hypothesis);

    // Simulating incremental addition of edges
    double avg_time = 0.0;
    cout << "Starting simulation of incremental dataset -> Displaying relative status : " << endl;
    for ( size_t candidate_id = 0 ; candidate_id < tot_hypothesis ; ++candidate_id )
    {
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        struct Info stat;

        // Selecting in deterministic way the candidate
        auto loop_candidate = gt_loops[candidate_id].second;

        // Start & End Idx of Simple Loop
        int id1 = loop_candidate->vertices()[0]->id();
        int id2 = loop_candidate->vertices()[1]->id();
        pair<int, int> candidate(id1, id2);
        vector<int> vcl_id = getClusterId(clusters, candidate);
        stat.gt = gt_loops[candidate_id].first;

        // No intersection found, creation of new cluster
        if ( vcl_id.size() == 0 )
        {
            // Adding the reliable voters
            OptimizableGraph::EdgeSet eset, eset_loop;
            for ( size_t j = id1 ; j < id2; eset.insert(odom_edges[j++]) );
            eset.insert(loop_candidate);

            // Generate map hypothesis closed-loop subproblem
            store(open_loop_problem);
            fixComplementary(open_loop_problem, id1, id2);
            stat.n_variables = id2 - id1;
            stat.n_loops = 0;

            if (!isAgreeingWithCurrentState(open_loop_problem, eset, cfg.fast_reject_th, cfg.fast_reject_iter_base, stat.opt_iterations))
            {
                outliers.push_back(loop_candidate);
                bucket[candidate_id] = 0;
                stat.estimate = false;
                restore(open_loop_problem);
            }
            else
            {
                stat.estimate = true;
                discard(open_loop_problem);
                clusters.push_back(candidate);
                edgesxcl.push_back(eset);
                eset_loop.insert(loop_candidate);
                edgesxcl_only_loops.push_back(eset_loop);
                max_consensus_set.push_back(candidate);
                propagateCurrentGuess(open_loop_problem, id2, odom_edges);
            }
        }
        else
        {
            // If the loop connects more clusters, then use the measurments from those clusters
            OptimizableGraph::EdgeSet eset, eset_loops, test;

            // Create a fake vertex for sorting the clusters & fusing cluster edges
            vector<pair<int, int>> sorted_intersect_clust;
            int cluster_edges = 0; int only_loops = 0;
            int connect_edges = 0; int extra_edges = 0;
            for ( size_t idx = 0; idx < vcl_id.size(); ++idx )
            {
                sorted_intersect_clust.push_back(clusters[vcl_id[idx]]);
                eset.insert(edgesxcl[vcl_id[idx]].begin(), edgesxcl[vcl_id[idx]].end());
                eset_loops.insert(edgesxcl_only_loops[vcl_id[idx]].begin(), edgesxcl_only_loops[vcl_id[idx]].end());
            }
            sort(sorted_intersect_clust.begin(), sorted_intersect_clust.end(), cmpSecond);
            cluster_edges = eset.size();
            only_loops = eset_loops.size();
            int only_odom = cluster_edges - only_loops;

            // Adding odom edges between connected clusters
            for ( size_t idx = 0; idx < sorted_intersect_clust.size() - 1; ++idx )
                for ( size_t j = sorted_intersect_clust[idx].second ; j < sorted_intersect_clust[idx + 1].first; eset.insert(odom_edges[j++]) );
            
            connect_edges = eset.size() - cluster_edges;
            int new_edges = eset.size();

            // New cluster beggining/ending
            int new_cl_start = sorted_intersect_clust[0].first;
            int new_cl_end = sorted_intersect_clust[sorted_intersect_clust.size() - 1].second;

            // Adding extra odom edges to the left and to the right
            for ( size_t idx = candidate.first ; idx < new_cl_start; eset.insert(odom_edges[idx++]) );
            for ( size_t idx = new_cl_end ; idx < candidate.second ; eset.insert(odom_edges[idx++]) );

            // Compute new cluster extremes
            new_cl_start = new_cl_start < candidate.first ? new_cl_start : candidate.first;
            new_cl_end = new_cl_end > candidate.second ? new_cl_end : candidate.second;
            extra_edges = eset.size() - new_edges;

            // Finally! Adding new_candidate edge!!
            eset.insert(loop_candidate);
            if ( vcl_id.size() > 1 )
                propagateGuess(open_loop_problem, 0, odom_edges.size(), odom_edges);

            // Generate map hypothesis closed-loop subproblem
            store(open_loop_problem); 
            stat.n_variables = new_cl_end - new_cl_start;
            stat.n_loops = only_loops; 
            fixComplementary(open_loop_problem, new_cl_start, new_cl_end);
            eset_loops.insert(loop_candidate);

            if (!isAgreeingWithCurrentState(open_loop_problem, eset, 
                                            cfg.slow_reject_th, 
                                            cfg.slow_reject_iter_base, 
                                            stat.opt_iterations))
            {
                restore(open_loop_problem);
                bucket[candidate_id] = 0;
                stat.estimate = false;
            }
            else
            {
                stat.estimate = true;
                discard(open_loop_problem);            
                max_consensus_set.push_back(candidate);
                propagateCurrentGuess(open_loop_problem, new_cl_end, odom_edges);

                if ( vcl_id.size() > 1 )
                {
                    vector<pair<int, int>> clusters_tmp;
                    vector<OptimizableGraph::EdgeSet> edgesxcl_tmp;
                    vector<OptimizableGraph::EdgeSet> edgesxcl_only_loops_tmp;

                    for ( size_t cidx = 0; cidx < clusters.size(); ++cidx )
                    {
                        // Check if out of the extremes of the cluster
                        bool out_right = clusters[cidx].second <= new_cl_start ? true : false;
                        bool out_left = clusters[cidx].first >= new_cl_end ? true : false;

                        if ( !(out_right || out_left) ) continue;

                        // Adding if independent
                        clusters_tmp.push_back(clusters[cidx]);
                        edgesxcl_tmp.push_back(edgesxcl[cidx]);
                        edgesxcl_only_loops_tmp.push_back(edgesxcl_only_loops[cidx]);
                    }
                    clusters = clusters_tmp;
                    edgesxcl = edgesxcl_tmp;
                    edgesxcl_only_loops = edgesxcl_only_loops_tmp;

                    clusters.push_back(make_pair(new_cl_start, new_cl_end));
                    edgesxcl.push_back(eset);
                    edgesxcl_only_loops.push_back(eset_loops);
                }
                else 
                {
                    clusters[vcl_id[0]].first = new_cl_start;
                    clusters[vcl_id[0]].second = new_cl_end;
                    edgesxcl[vcl_id[0]] = eset;
                    edgesxcl_only_loops[vcl_id[0]] = eset_loops;
                }
                
            }

        }

        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        chrono::microseconds delta_time = chrono::duration_cast<chrono::microseconds>(end - begin);
        stat.opt_time = delta_time.count() / 1000000.0;
        avg_time += stat.opt_time;
        statistics[candidate_id] = stat;

        printProgress((double)(candidate_id + 1) / (double)tot_hypothesis);
    }
    cout << "\nCompleted!" << endl;

    propagateGuess(open_loop_problem, 0, odom_edges.size(), odom_edges);
    OptimizableGraph::EdgeSet eset_gl;
    for ( size_t i = 0 ; i < odom_edges.size(); eset_gl.insert(odom_edges[i++]) )
        odom_edges[i]->setInformation(odom_edges[i]->information() * 0.1);
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

    for ( size_t counter = 0 ; counter < clusters.size() ; ++counter )    
        cout << "L("<< counter << ") = [" << clusters[counter].first << ", " << clusters[counter].second << "]\n";
    cout << "Avg Time x test = " << avg_time / tot_hypothesis << " [s]\n";
    cout << "Precision = " << precision << endl;        
    cout << "Recall = " << recall << endl;        
    std::cout << "Size of MAX consistent set = " << max_consensus_set.size() << std::endl;

    generateDistributionFile(cfg.name, statistics);

    ofstream outfile;
    outfile.open(output_file.c_str());
    for (size_t it = 0; it < open_loop_problem.vertices().size(); ++it )
    {
        VertexSE2* v = dynamic_cast<VertexSE2*>(open_loop_problem.vertex(it));
        outfile << v->estimate()[0] << " " << v->estimate()[1] << " " << v->estimate()[2] << endl;
    }
    outfile.close();


    string out2 = output_file.substr(0, output_file.size() - 3) + "PR";
    outfile.open(out2.c_str());
    outfile << precision << " " << recall << endl;
    outfile << dt << endl;
    outfile.close();

    return;
}
