#include "ipc/consensus.hpp"
#include "switchableConstraints/include/edge_switchPrior.hpp"

using namespace std;
using namespace g2o;
using namespace Eigen;

// Implementing naive recovery with SC
template <class EDGE, class VERTEX> 
void IPC<EDGE, VERTEX>::recoverySC(std::vector<int>& in_out, 
                                   const std::vector<EDGE*>& loops)
{   
    // Assigning different weights based on the results of IPC
    double weight_in = 0.9;
    double weight_out = 0.2;
    double information_edge = 10.0;
    int id_swv = _problem->vertices().size() + 300;

    OptimizableGraph::EdgeSet eset;
    vector<VertexSwitchLinear*> recovery_vertices;
    for (int candidate_id = 0; candidate_id < loops.size(); ++candidate_id)
    {
        EDGE* loop_candidate = loops[candidate_id];
        double weight = in_out[candidate_id] == 1 ? weight_in : weight_out;
        
        // Create a switchable vertex
        VertexSwitchLinear* sw = new VertexSwitchLinear();
        sw->setEstimate(weight);
        sw->setFixed(false);
        sw->setId(id_swv);
        recovery_vertices.push_back(sw);

        // Prior edge based on the results of IPC
        EdgeSwitchPrior* sw_prior = new EdgeSwitchPrior();
        sw_prior->setVertex(0, sw);
        Eigen::Matrix<double, 1, 1> information_mat = Eigen::Matrix<double, 1, 1>::Identity();
        sw_prior->setMeasurement(weight);
        sw_prior->setInformation(information_mat * information_edge);    

        // Switchable edge
        auto sw_constraint = getSwitchableEdge(loop_candidate, sw, weight);

        // Add everything to the problem
        _problem->addVertex(sw);
        _problem->addEdge(sw_prior);
        _problem->addEdge(sw_constraint);
    
        // Create a set of edges to be optimized
        eset.insert(sw_constraint);
        eset.insert(sw_prior);

        ++id_swv;
    }

    // Add the odometry edges to the optimization set
    for ( size_t o_id = 0 ; o_id < _odom_edges.size(); eset.insert(_odom_edges[o_id++]) );

    // Optimize the problem
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    _problem->initializeOptimization(eset);
    _problem->optimize(1000);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    chrono::microseconds delta_time = chrono::duration_cast<chrono::microseconds>(end - begin);
    double dt = delta_time.count() / 1000000.0;

    // Updating values 
    float agree = 0.0;
    for (int candidate_id = 0; candidate_id < loops.size(); ++candidate_id)
    {
        VertexSwitchLinear* sw = recovery_vertices[candidate_id];
        int bucket_val = sw->estimate() > 0.5 ? 1 : 0;

        agree += bucket_val == in_out[candidate_id] ? 1 : 0;
        in_out[candidate_id] = bucket_val;
    }
    agree /= loops.size();
    cout << "Agreement after recovery: " << agree * 100.0 << endl;
    cout << "Disagreement after recovery: " << (1.0 - agree) * 100.0 << endl;
    cout << "Time for recovery: " << dt << " seconds" << endl;

    return;
}

template void IPC<EdgeSE3, VertexSE3>::recoverySC(std::vector<int>& in_out, const std::vector<EdgeSE3*>& loops);
template void IPC<EdgeSE2, VertexSE2>::recoverySC(std::vector<int>& in_out, const std::vector<EdgeSE2*>& loops);

/*----------------------------------------------------------------*/

template <class EDGE, class VERTEX>
EdgeSE3Switchable* IPC<EDGE, VERTEX>::getSwitchableEdge(EdgeSE3* candidate, VertexSwitchLinear* sw, double weight)
{
    EdgeSE3Switchable* edge = new EdgeSE3Switchable();
    edge->setVertex(0, candidate->vertex(0));
    edge->setVertex(1, candidate->vertex(1));
    edge->setVertex(2, sw);
    Eigen::Isometry3d meas = candidate->measurement(); 
    edge->setMeasurement(SE3Quat(meas.rotation(), meas.translation()));
    edge->setInformation(candidate->information());

    return edge;
}

template <class EDGE, class VERTEX>
EdgeSE2Switchable* IPC<EDGE, VERTEX>::getSwitchableEdge(EdgeSE2* candidate, VertexSwitchLinear* sw, double weight)
{
    EdgeSE2Switchable* edge = new EdgeSE2Switchable();
    edge->setVertex(0, candidate->vertex(0));
    edge->setVertex(1, candidate->vertex(1));
    edge->setVertex(2, sw);
    edge->setMeasurement(candidate->measurement());
    edge->setInformation(candidate->information());

    return edge;
}


template EdgeSE3Switchable* IPC<EdgeSE3, VertexSE3>::getSwitchableEdge(EdgeSE3* candidate, VertexSwitchLinear* sw, double weight);
template EdgeSE2Switchable* IPC<EdgeSE2, VertexSE2>::getSwitchableEdge(EdgeSE2* candidate, VertexSwitchLinear* sw, double weight);

/*----------------------------------------------------------------*/