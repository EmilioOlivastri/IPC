#pragma once

#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/core/base_multi_edge.h"
#include "g2o/core/hyper_graph_action.h"

class EdgeSE2Switchable : public g2o::BaseMultiEdge<3, g2o::SE2>
{
  public:
    EdgeSE2Switchable();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    void computeError();
    void linearizeOplus();

};