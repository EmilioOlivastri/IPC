#pragma once

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/base_multi_edge.h"
#include "g2o/core/hyper_graph_action.h"

class EdgeSE3Switchable : public g2o::BaseMultiEdge<6, g2o::SE3Quat>
{
  public:
    EdgeSE3Switchable();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    void computeError();
    void linearizeOplus();

};