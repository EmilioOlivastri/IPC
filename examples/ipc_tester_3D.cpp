#include "ipc/consensus.hpp"
#include "ipc/simulation.hpp"

using namespace std;
using namespace g2o;
using namespace Eigen;

G2O_USE_TYPE_GROUP(slam3d);
G2O_USE_OPTIMIZATION_LIBRARY(eigen);

void correctedInformationMatrices(g2o::SparseOptimizer& optimizer);

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
  correctedInformationMatrices(optimizer);

  vector<EdgeSE3*> loops, odom_edges;
  splitProblemConstraints<EdgeSE3>(optimizer, odom_edges, loops);
  simulating_incremental_data<Isometry3d, EdgeSE3, VertexSE3>(cfg, optimizer, loops);

  return 0;
}


void correctedInformationMatrices(g2o::SparseOptimizer& optimizer)
{
  for ( auto it_e = optimizer.edges().begin(); it_e != optimizer.edges().end(); ++it_e )
  {
      auto edge = dynamic_cast<EdgeSE3*>(*it_e);
      if ( edge == nullptr ) continue;

      Eigen::Matrix<double, 6, 6> info = Eigen::Matrix<double, 6, 6>::Zero();

      // Correcting translation part
      info(0, 0) = 800.0;  // TUM: 300.0   | KITTI_05: 1000.0 | KITTI_00: 800.0
      info(1, 1) = 650.0;  // TUM: 260.0   | KITTI_05: 800.0  | KITTI_00: 650.0
      info(2, 2) = 5000.0; // TUM: 2500.0  | KITTI_05: 1000.0 | KITTI_00: 5000.0
      // Correcting rotation part
      info(3, 3) = 5000.0; // TUM: 2500.0 | KITTI_05: 1000.0 | KITTI_00: 5000.0
      info(4, 4) = 5000.0; // TUM: 2500.0 | KITTI_05: 1000.0 | KITTI_00: 5000.0
      info(5, 5) = 700.0;  // TUM: 200.0  | KITTI_05: 900.0  | KITTI_00: 700.0
      edge->setInformation(info);
  }  

  return;
}