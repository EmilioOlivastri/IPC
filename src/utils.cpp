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

// Store the state of the graph
bool store(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VertexSE2*>(vIt->second)->push();

    return true;
}

bool restore(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VertexSE2*>(vIt->second)->pop();

    store(opt);

    return true;
}

bool discard(SparseOptimizer& opt)
{
    for (auto vIt = opt.vertices().begin(); vIt != opt.vertices().end(); ++vIt)
        static_cast<VertexSE2*>(vIt->second)->discardTop();

    return true;
}

void odometryInitialization(SparseOptimizer& optimizer)
{
    // Iterating trough the edges
    for ( auto it_e = optimizer.edges().begin(); it_e != optimizer.edges().end(); ++it_e )
    {
        auto edge_odom = dynamic_cast<EdgeSE2*>(*it_e);
        if ( edge_odom != nullptr )
        {
            auto v_fn = dynamic_cast<VertexSE2*>(edge_odom->vertices()[1]);
            auto v_st = dynamic_cast<VertexSE2*>(edge_odom->vertices()[0]);

            int id_st = v_st->id();
            int id_fn = v_fn->id();

            auto est = id_fn - id_st == 1 ? v_st->estimate() * edge_odom->measurement(): v_fn->estimate();
            v_fn->setEstimate(est);
        }
    }

    return;
}

void setProblem(const string& problem_file, 
                SparseOptimizer& optimizer,
                vector<SE2>& init_poses,
                vector<VertexSE2*>& v_poses)
{
    optimizer.setVerbose(false);
    
    // allocate the solver
    g2o::OptimizationAlgorithmProperty solverProperty;
    optimizer.setAlgorithm(g2o::OptimizationAlgorithmFactory::instance()->construct("dl_var", solverProperty));
    
    // Loading the g2o file
    ifstream ifs(problem_file.c_str());
    if (!ifs) 
    {
        cerr << "unable to open " << problem_file << endl;
        return;
    }
    optimizer.load(ifs);
    odometryInitialization(optimizer);
    init_poses.resize(optimizer.vertices().size());
    v_poses.resize(optimizer.vertices().size());
    for ( auto it = optimizer.vertices().begin() ; it != optimizer.vertices().end() ; ++it )
    {
        auto v1 = dynamic_cast<VertexSE2*>(it->second);
        v_poses[v1->id()] = v1;
        init_poses[v1->id()]= v1->estimate();
    }

    return;
}

void getProblemNOLOOPS(SparseOptimizer& optimizer,
                       vector<EdgeSE2*>& loops)
{
    // Removing all the edges of the problem
    vector<EdgeSE2*> toremove;
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EdgeSE2*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) > 1 )
        {
            EdgeSE2* edge_cpy = new EdgeSE2;
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


void splitProblemConstraints(SparseOptimizer& optimizer,
                             vector<EdgeSE2*>& odom,
                             vector<EdgeSE2*>& loops)
{
    // Dividing all the edges of the problem
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EdgeSE2*>(*it);
        if ( edge == nullptr) continue;
        
        int diff = abs(edge->vertices()[1]->id() - edge->vertices()[0]->id());

        diff > 1 ? loops.push_back(edge) : odom.push_back(edge);
    }

    return;
}

// Controllato e giusto
void getProblemLoops(SparseOptimizer& optimizer,
                     vector<EdgeSE2*>& loops)
{
    // Getting the loops edges
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EdgeSE2*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) > 1 )
            loops.push_back(edge);
    }

    return;
}

// Controllato e giusto
void getProblemOdom(SparseOptimizer& optimizer,
                    vector<EdgeSE2*>& odom)
{
    // Getting the loops edges
    for ( auto it = optimizer.edges().begin() ; it != optimizer.edges().end() ; ++it )
    {
        auto edge = dynamic_cast<EdgeSE2*>(*it);
        if ( edge != nullptr && abs(edge->vertices()[1]->id() - edge->vertices()[0]->id()) == 1 )
            odom.push_back(edge);
    }
}

void store(const std::string& name, const std::vector<g2o::VertexSE2*>& v_poses)
{
    cout << "Saving final poses to file: " << name << endl;
    ofstream poses_file(name.c_str());

    // Iterating the vector of vertices to extract the map
    for ( int i = 0 ; i < v_poses.size() ; ++i )
    {
        SE2 opt_pose = v_poses[i]->estimate(); 
        poses_file << opt_pose[0] << " " << opt_pose[1] << " " << opt_pose[2] << endl;
    }

    poses_file.close();

    return;
}

void readSolutionFile(vector<Eigen::Isometry2d>& poses, const string& path)
{
    const int LINESIZE = 81920;

    std::ifstream in_data(path.c_str());
    if (!in_data )
        throw std::invalid_argument("Cannot find file : " + path);

    // Keep going until the file has been read
    double x, y, yaw;
    while ( !in_data.eof() )
    {
        in_data >> x >> y >> yaw;

        Eigen::Isometry2d p = Eigen::Isometry2d::Identity();
        p.translation() = Eigen::Vector2d(x, y);
        p.linear() = Eigen::Rotation2D<double>(yaw).toRotationMatrix();

        poses.push_back(p);

        in_data.ignore(LINESIZE, '\n');
    }

    return;
}

void readConfig(const std::string& cfg_filepath, Config& out_cfg)
{
    const YAML::Node config = YAML::LoadFile(cfg_filepath);

    // Filter parameters
    out_cfg.name = config["name"].as<std::string>();
    out_cfg.dataset = config["dataset"].as<std::string>();
    out_cfg.ground_truth = config["ground_truth"].as<std::string>();
    out_cfg.output = config["output"].as<std::string>();
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

bool cmpEdgesID(EdgeSE2* e1, EdgeSE2* e2)
{
    return (e1->vertices()[1]->id() < e2->vertices()[1]->id());
}

bool cmpTime(pair<int, EdgeSE2*> p1, pair<int, EdgeSE2*> p2)
{
    int id1_v1 = p1.second->vertices()[1]->id();
    int id1_v2 = p1.second->vertices()[0]->id();
    int id1_max = id1_v1 > id1_v2 ? id1_v1 : id1_v2;

    int id2_v1 = p2.second->vertices()[1]->id();
    int id2_v2 = p2.second->vertices()[0]->id();
    int id2_max = id2_v1 > id2_v2 ? id2_v1 : id2_v2;

    return (id1_max < id2_max);
}