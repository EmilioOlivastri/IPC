#pragma once 

#include <iostream>
#include <yaml-cpp/yaml.h>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/command_args.h"
#include "g2o/types/slam2d/types_slam2d.h"
#include "g2o/types/slam3d/types_slam3d.h"

struct Config
{
  std::string name;
  std::string dataset;
  std::string ground_truth;
  std::string output;
  bool visualize;
  int canonic_inliers;
  double fast_reject_th;
  int fast_reject_iter_base;
  double slow_reject_th;
  int slow_reject_iter_base;
};


// Using edges to initialize graph
template <class EDGE, class VERTEX>
void odometryInitialization(g2o::SparseOptimizer& optimizer);

// Sets the whole optimization problem
template <class T, class EDGE, class VERTEX>
void setProblem(const std::string& problem_file, 
                g2o::SparseOptimizer& optimizer,
                std::vector<T>& init_poses,
                std::vector<VERTEX*>& v_poses);

// Gets the problem with no loops to then test them singularly
template <class EDGE>
void getProblemNOLOOPS(g2o::SparseOptimizer& optimizer,
                       std::vector<EDGE*>& loops);

// Get the loops from the problem without explicitly eliminating them
template <class EDGE>
void getProblemLoops(g2o::SparseOptimizer& optimizer,
                     std::vector<EDGE*>& loops);

// Get only the odometry edges
template <class EDGE>
void splitProblemConstraints(g2o::SparseOptimizer& optimizer,
                             std::vector<EDGE*>& odom,
                             std::vector<EDGE*>& loops);

// Get only the odometry edges
template <class EDGE>
void getProblemOdom(g2o::SparseOptimizer& optimizer,
                    std::vector<EDGE*>& odom);

// Store the state of the graph
template <class VERTEX>
bool store(g2o::SparseOptimizer& opt);

// Restore the state of the graph
template <class VERTEX>
bool restore(g2o::SparseOptimizer& opt);

// Discards the last estimate
template <class VERTEX>
bool discard(g2o::SparseOptimizer& opt);

void writeVertex(std::ofstream& out_data, g2o::VertexSE2* v);
void writeVertex(std::ofstream& out_data, g2o::VertexSE3* v);

template <class T>
void readSolutionFile(std::vector<T>& poses, const std::string& path);
void readLine(std::ifstream& in_data, Eigen::Isometry2d& pose); 
void readLine(std::ifstream& in_data, Eigen::Isometry3d& pose);
void readConfig(const std::string& cfg_filepath, Config& out_cfg);
void printProgress(double percentage);

/*
Compare functions for the sorting/comparing operations from STL
*/
bool cmpFirst(std::pair<int, int> p1, std::pair<int, int> p2);
bool cmpSecond(std::pair<int, int> p1, std::pair<int, int> p2);
bool cmpEdgesID(g2o::OptimizableGraph::Edge* e1, 
                g2o::OptimizableGraph::Edge* e2);
bool cmpTime(std::pair<int, g2o::OptimizableGraph::Edge*> p1, 
             std::pair<int, g2o::OptimizableGraph::Edge*> p2);