#include "ipc/consensus.hpp"

using namespace std;
using namespace g2o;

G2O_USE_TYPE_GROUP(slam2d);
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
  vector<SE2> init_poses;
  vector<VertexSE2*> v_poses;

  // create the optimizer to load the data and carry out the optimization
  SparseOptimizer optimizer;
  setProblem(cfg.dataset, optimizer, init_poses, v_poses);

  vector<EdgeSE2*> loops, odom_edges;
  splitProblemConstraints(optimizer, odom_edges, loops);

  simulating_incremental_data(cfg, optimizer, loops);

  return 0;
}