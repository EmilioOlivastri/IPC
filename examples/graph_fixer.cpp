#include "ipc/consensus.hpp"

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

  cout << "Tot Edges = " << optimizer.edges().size() << endl;
  cout << "Tot Vertices = " << optimizer.vertices().size() << endl;

  vector<EdgeSE3*> loops, odom_edges;
  splitProblemConstraints<EdgeSE3>(optimizer, odom_edges, loops);

  for (size_t idx = 0; idx < loops.size(); ++idx)
  {
    // If v1->id() < v2->id() then it is in correct order
    if ( loops[idx]->vertices()[0]->id() < loops[idx]->vertices()[1]->id() )
        continue;

    EdgeSE3* new_loop = new EdgeSE3;
    new_loop->vertices()[0] = loops[idx]->vertices()[1];
    new_loop->vertices()[1] = loops[idx]->vertices()[0];
    new_loop->setInformation(loops[idx]->information());
    new_loop->setMeasurement(loops[idx]->measurement().inverse());

    optimizer.addEdge(new_loop);
    optimizer.removeEdge(loops[idx]);
  }
  
  optimizer.save("graph.g2o");

  optimizer.clear();


  return 0;
}