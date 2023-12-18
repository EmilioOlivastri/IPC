#include "voting_consensus/utils.hpp"
#include "voting_consensus/consensus.hpp"
#include "voting_consensus/consensus_utils.hpp"

using namespace std;
using namespace g2o;

// we use the 2D and 3D SLAM types here
G2O_USE_TYPE_GROUP(slam2d);
G2O_USE_TYPE_GROUP(slam3d);
G2O_USE_OPTIMIZATION_LIBRARY(eigen);

int main(int argc, char** argv) 
{
  // Command line parsing
  string cfgFilename;
  string inputFilename;
  string outputFilename;
  CommandArgs arg;
  arg.param("c", cfgFilename, "",
            "path to cfg file");
  arg.param("o", outputFilename, "", "output filename");
  arg.paramLeftOver("graph-input", inputFilename, "",
                    "graph file which will be processed");
  arg.parseArgs(argc, argv);

  Config cfg;
  readConfig(cfgFilename, cfg);

  // Storing initial guess and vertices of optimization
  vector<SE2> init_poses;
  vector<VertexSE2*> v_poses;

  // create the optimizer to load the data and carry out the optimization
  SparseOptimizer optimizer;
  setProblem(inputFilename, optimizer, init_poses, v_poses);

  vector<EdgeSE2*> loops, odom_edges;
  splitProblemConstraints(optimizer, odom_edges, loops);

  vector<EdgeSE2*> inliers, outliers;
  ipc(cfg, optimizer, loops, inliers, outliers, outputFilename);

  return 0;
}