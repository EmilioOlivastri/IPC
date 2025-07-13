#include "ipc/consensus.hpp"
#include "ipc/simulation.hpp"

using namespace std;
using namespace g2o;
using namespace Eigen;

G2O_USE_TYPE_GROUP(slam3d);
G2O_USE_OPTIMIZATION_LIBRARY(eigen);

int main(int argc, char** argv) 
{
  // Command line parsing
  string cfgFilename;
  CommandArgs arg;
  arg.param("c", cfgFilename, "",
            "path to cfg file");
  arg.parseArgs(argc, argv);

  Config cfg;
  readConfig(cfgFilename, cfg);

  // Storing initial guess and vertices of optimization
  vector<Isometry3d> init_poses;
  vector<VertexSE3*> v_poses;

  // create the optimizer to load the data and carry out the optimization
  SparseOptimizer optimizer;
  setProblem<Isometry3d, EdgeSE3, VertexSE3>(cfg.dataset, optimizer, init_poses, v_poses);

  vector<EdgeSE3*> loops, odom_edges;
  splitProblemConstraints<EdgeSE3>(optimizer, odom_edges, loops);
  simulating_incremental_data<Isometry3d, EdgeSE3, VertexSE3>(cfg, optimizer, loops);

  return 0;
}