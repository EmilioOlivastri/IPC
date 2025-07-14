#include "ipc/consensus.hpp"
#include "ipc/simulation.hpp"

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

        //printProgress((double)(candidate_id + 1) / (double)tot_hypothesis);
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

    //const vector<pair<int, int>> clusters = ipc.getClusters();
    const vector<EDGE*> max_consensus_set = ipc.getMaxConsensusSet();
    cout << "Size of MAX consistent set = " << max_consensus_set.size() << endl;

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