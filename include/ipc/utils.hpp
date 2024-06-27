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
void odometryInitialization(g2o::SparseOptimizer& optimizer);

// Sets the whole optimization problem
void setProblem(const std::string& problem_file, 
                g2o::SparseOptimizer& optimizer,
                std::vector<g2o::SE2>& init_poses,
                std::vector<g2o::VertexSE2*>& v_poses);

// Gets the problem with no loops to then test them singularly
void getProblemNOLOOPS(g2o::SparseOptimizer& optimizer,
                       std::vector<g2o::EdgeSE2*>& loops);

// Get the loops from the problem without explicitly eliminating them
void getProblemLoops(g2o::SparseOptimizer& optimizer,
                     std::vector<g2o::EdgeSE2*>& loops);

// Get only the odometry edges
void splitProblemConstraints(g2o::SparseOptimizer& optimizer,
                             std::vector<g2o::EdgeSE2*>& odom,
                             std::vector<g2o::EdgeSE2*>& loops);

// Get only the odometry edges
void getProblemOdom(g2o::SparseOptimizer& optimizer,
                    std::vector<g2o::EdgeSE2*>& odom);

// Store the state of the graph
bool store(g2o::SparseOptimizer& opt);

// Restore the state of the graph
bool restore(g2o::SparseOptimizer& opt);

// Discards the last estimate
bool discard(g2o::SparseOptimizer& opt);
void readSolutionFile(std::vector<Eigen::Isometry2d>& poses, const std::string& path);
void readConfig(const std::string& cfg_filepath, Config& out_cfg);
void printProgress(double percentage);

/*
Compare functions for the sorting/comparing operations from STL
*/
bool cmpFirst(std::pair<int, int> p1, std::pair<int, int> p2);
bool cmpSecond(std::pair<int, int> p1, std::pair<int, int> p2);
bool cmpEdgesID(g2o::EdgeSE2* e1, g2o::EdgeSE2* e2);
bool cmpTime(std::pair<int, g2o::EdgeSE2*> p1, std::pair<int, g2o::EdgeSE2*> p2);