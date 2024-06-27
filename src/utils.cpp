#include "ipc/utils.hpp"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/solvers/pcg/linear_solver_pcg.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/core/optimization_algorithm_dogleg.h"

#include <string>
#include <iostream>
#include <fstream>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace std;
using namespace g2o;
using namespace Eigen;

// Store the state of the graph
template <class VERTEX>
bool store(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VERTEX*>(vIt->second)->push();

    return true;
}

template bool store<VertexSE2>(SparseOptimizer& opt);
template bool store<VertexSE3>(SparseOptimizer& opt);

/*-------------------------------------------------------------------*/

template <class VERTEX>
bool restore(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VERTEX*>(vIt->second)->pop();

    store<VERTEX>(opt);

    return true;
}

template bool restore<VertexSE2>(SparseOptimizer& opt);
template bool restore<VertexSE3>(SparseOptimizer& opt);

/*-------------------------------------------------------------------*/

template <class VERTEX>
bool discard(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VERTEX*>(vIt->second)->discardTop();

    return true;
}

template bool discard<VertexSE2>(SparseOptimizer& opt);
template bool discard<VertexSE3>(SparseOptimizer& opt);

/*-------------------------------------------------------------------*/

template <class EDGE, class VERTEX>
void odometryInitialization(SparseOptimizer& optimizer)
{
    // Iterating trough the edges
    for ( auto it_e = optimizer.edges().begin(); it_e != optimizer.edges().end(); ++it_e )
    {
        auto edge_odom = dynamic_cast<EDGE*>(*it_e);
        if ( edge_odom != nullptr )
        {
            auto v_fn = dynamic_cast<VERTEX*>(edge_odom->vertices()[1]);
            auto v_st = dynamic_cast<VERTEX*>(edge_odom->vertices()[0]);

            int id_st = v_st->id();
            int id_fn = v_fn->id();

            auto est = id_fn - id_st == 1 ? v_st->estimate() * edge_odom->measurement(): v_fn->estimate();
            v_fn->setEstimate(est);
        }
    }

    return;
}

template void odometryInitialization<EdgeSE2, VertexSE2>(SparseOptimizer& optimizer);
template void odometryInitialization<EdgeSE3, VertexSE3>(SparseOptimizer& optimizer);

/*-------------------------------------------------------------------*/

template <class T, class EDGE, class VERTEX>
void setProblem(const string& problem_file, 
                SparseOptimizer& optimizer,
                vector<T>& init_poses,
                vector<VERTEX*>& v_poses)
{
    optimizer.setVerbose(false);
    
    // allocate the solver
    OptimizationAlgorithmProperty solverProperty;
    optimizer.setAlgorithm(OptimizationAlgorithmFactory::instance()->construct("dl_var", solverProperty));
    
    // Loading the g2o file
    ifstream ifs(problem_file.c_str());
    if (!ifs) 
    {
        cerr << "unable to open " << problem_file << endl;
        return;
    }
    optimizer.load(ifs);
    odometryInitialization<EDGE, VERTEX>(optimizer);
    init_poses.resize(optimizer.vertices().size());
    v_poses.resize(optimizer.vertices().size());
    for ( auto it = optimizer.vertices().begin() ; it != optimizer.vertices().end() ; ++it )
    {
        auto v1 = dynamic_cast<VERTEX*>(it->second);
        v_poses[v1->id()] = v1;
        init_poses[v1->id()]= v1->estimate();
    }

    return;
}

template void setProblem<SE2, EdgeSE2, VertexSE2>(const string& problem_file, 
                                                  SparseOptimizer& optimizer,
                                                  vector<SE2>& init_poses,
                                                  vector<VertexSE2*>& v_poses);
template void setProblem<Isometry3d, EdgeSE3, VertexSE3>(const string& problem_file, 
                                                         SparseOptimizer& optimizer,
                                                         vector<Isometry3d>& init_poses,
                                                         vector<VertexSE3*>& v_poses);

/*-------------------------------------------------------------------*/

template <class EDGE>
void getProblemNOLOOPS(SparseOptimizer& optimizer,
                       vector<EDGE*>& loops)
{
    // Removing all the edges of the problem
    vector<EDGE*> toremove;
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EDGE*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) > 1 )
        {
            EDGE* edge_cpy = new EDGE;
            edge_cpy->vertices()[0] = edge->vertices()[0];
            edge_cpy->vertices()[1] = edge->vertices()[1];
            edge_cpy->setMeasurement(edge->measurement());
            edge_cpy->setInformation(edge->information());
            loops.push_back(edge_cpy);
            toremove.push_back(edge);
        }
    }

    // Needs to remove using second loop, otherwise it breake the loop above
    for ( size_t i = 0 ; i < toremove.size() ; ++i )
        optimizer.removeEdge(toremove[i]);

    return;
}

template void getProblemNOLOOPS<EdgeSE2>(SparseOptimizer& optimizer, vector<EdgeSE2*>& loops);
template void getProblemNOLOOPS<EdgeSE3>(SparseOptimizer& optimizer, vector<EdgeSE3*>& loops);

/*-------------------------------------------------------------------*/

template <class EDGE>
void splitProblemConstraints(SparseOptimizer& optimizer,
                             vector<EDGE*>& odom,
                             vector<EDGE*>& loops)
{
    // Dividing all the edges of the problem
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EDGE*>(*it);
        if ( edge == nullptr) continue;
        
        int diff = abs(edge->vertices()[1]->id() - edge->vertices()[0]->id());

        diff > 1 ? loops.push_back(edge) : odom.push_back(edge);
    }

    return;
}

template void splitProblemConstraints<EdgeSE2>(SparseOptimizer& optimizer, vector<EdgeSE2*>& odom, vector<EdgeSE2*>& loops);
template void splitProblemConstraints<EdgeSE3>(SparseOptimizer& optimizer, vector<EdgeSE3*>& odom, vector<EdgeSE3*>& loops);

/*-------------------------------------------------------------------*/

// Controllato e giusto
template <class EDGE>
void getProblemLoops(SparseOptimizer& optimizer,
                     vector<EDGE*>& loops)
{
    // Getting the loops edges
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EDGE*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) > 1 )
            loops.push_back(edge);
    }

    return;
}

template void getProblemLoops<EdgeSE2>(SparseOptimizer& optimizer, vector<EdgeSE2*>& loops);
template void getProblemLoops<EdgeSE3>(SparseOptimizer& optimizer, vector<EdgeSE3*>& loops);

/*-------------------------------------------------------------------*/

// Controllato e giusto
template <class EDGE>
void getProblemOdom(SparseOptimizer& optimizer,
                    vector<EDGE*>& odom)
{
    // Getting the loops edges
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EDGE*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) == 1 )
            odom.push_back(edge);
    }

    return;
}

template void getProblemOdom<EdgeSE2>(SparseOptimizer& optimizer, vector<EdgeSE2*>& odom);
template void getProblemOdom<EdgeSE3>(SparseOptimizer& optimizer, vector<EdgeSE3*>& odom);


/*-------------------------------------------------------------------*/

void writeVertex(ofstream& out_data, VertexSE2* v)
{
    out_data << v->estimate()[0] << " " 
             << v->estimate()[1] << " " 
             << v->estimate()[2] << endl;
    return;
}


void writeVertex(ofstream& out_data, VertexSE3* v)
{
    Isometry3d pose = v->estimate();
    Quaterniond q(pose.linear());
    Vector3d t = pose.translation();

    out_data << t[0] << " " << t[1] << " " << t[2] << " "
             << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << endl;

    return;
}

/*-------------------------------------------------------------------*/

template <class T>
void readSolutionFile(vector<T>& poses, const string& path)
{
    const int LINESIZE = 81920;

    ifstream in_data(path.c_str());
    if (!in_data )
        throw invalid_argument("Cannot find file : " + path);

    // Keep going until the file has been read
    while ( !in_data.eof() )
    {
        T p = T::Identity();
        readLine(in_data, p);
        poses.push_back(p);

        in_data.ignore(LINESIZE, '\n');
    }

    return;
}

template void readSolutionFile<Isometry2d>(vector<Isometry2d>& poses, const string& path);
template void readSolutionFile<Isometry3d>(vector<Isometry3d>& poses, const string& path);

/*-------------------------------------------------------------------*/

void readLine(ifstream& in_data, Isometry2d& pose)
{
    double x, y, yaw;
    in_data >> x >> y >> yaw;

    pose = Isometry2d::Identity();
    pose.translation() = Vector2d(x, y);
    pose.linear() = Eigen::Rotation2D<double>(yaw).toRotationMatrix();

    return;
}

void readLine(ifstream& in_data, Isometry3d& pose)
{
    double x, y, z, qx, qy, qz, qw;
    in_data >> x >> y >> z >> qx >> qy >> qz >> qw;

    pose = Isometry3d::Identity();
    pose.translation() = Vector3d(x, y, z);
    pose.linear() = Quaterniond(qw, qx, qy, qz).toRotationMatrix();

    return;
}


/*-------------------------------------------------------------------*/

void readConfig(const string& cfg_filepath, Config& out_cfg)
{
    const YAML::Node config = YAML::LoadFile(cfg_filepath);

    // Filter parameters
    out_cfg.name = config["name"].as<string>();
    out_cfg.dataset = config["dataset"].as<string>();
    out_cfg.ground_truth = config["ground_truth"].as<string>();
    out_cfg.output = config["output"].as<string>();
    out_cfg.visualize = config["visualize"].as<int>() == 1 ? true : false;
    out_cfg.canonic_inliers = config["canonic_inliers"].as<int>();
    out_cfg.fast_reject_th = config["fast_reject_th"].as<double>();
    out_cfg.fast_reject_iter_base = config["fast_reject_iter_base"].as<int>();
    out_cfg.slow_reject_th = config["slow_reject_th"].as<double>();
    out_cfg.slow_reject_iter_base = config["slow_reject_iter_base"].as<int>();

    return;
}

void printProgress(double percentage) 
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

bool cmpFirst(pair<int, int> p1, pair<int, int> p2)
{
    return (p1.first < p2.first);
}


bool cmpSecond(pair<int, int> p1, pair<int, int> p2)
{
    return (p1.second < p2.second);
}

bool cmpEdgesID(OptimizableGraph::Edge* e1, OptimizableGraph::Edge* e2)
{
    return (e1->vertices()[1]->id() < e2->vertices()[1]->id());
}

bool cmpTime(pair<int, OptimizableGraph::Edge*> p1, pair<int, OptimizableGraph::Edge*> p2)
{
    int id1_v1 = p1.second->vertices()[1]->id();
    int id1_v2 = p1.second->vertices()[0]->id();
    int id1_max = id1_v1 > id1_v2 ? id1_v1 : id1_v2;

    int id2_v1 = p2.second->vertices()[1]->id();
    int id2_v2 = p2.second->vertices()[0]->id();
    int id2_max = id2_v1 > id2_v2 ? id2_v1 : id2_v2;

    return (id1_max < id2_max);
}