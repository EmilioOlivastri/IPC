#include "edge_se3Switchable.hpp"
#include "vertex_switchLinear.hpp"
#include "g2o/types/slam3d/se3quat.h"
#include "jacobian.hpp"

using namespace std;
using namespace Eigen;


// ================================================
EdgeSE3Switchable::EdgeSE3Switchable() : g2o::BaseMultiEdge<6, g2o::SE3Quat>()
{
  resize(3);
  _jacobianOplus[0].resize(6,6); 
  _jacobianOplus[1].resize(6,6);
  _jacobianOplus[2].resize(6,1);

}
// ================================================
bool EdgeSE3Switchable::read(std::istream& is)
  {
    return true;
  }
// ================================================
bool EdgeSE3Switchable::write(std::ostream& os) const
{
    return os.good();
}

// ================================================
void EdgeSE3Switchable::linearizeOplus()
{

    g2o::VertexSE3* from = static_cast<g2o::VertexSE3*>(_vertices[0]);
    g2o::VertexSE3* to = static_cast<g2o::VertexSE3*>(_vertices[1]);
    const VertexSwitchLinear* vSwitch = static_cast<const VertexSwitchLinear*>(_vertices[2]);

    Matrix3d izR        = (measurement().inverse()).rotation().toRotationMatrix();
    const Vector3d& izt = (measurement().inverse()).translation();

    Eigen::Isometry3d iXiXj_iso = from->estimate().inverse() * to->estimate();
    g2o::SE3Quat iXiXj         = g2o::SE3Quat(iXiXj_iso.rotation(), iXiXj_iso.translation());
    Matrix3d iRiRj        = iXiXj.rotation().toRotationMatrix();
    const Vector3d& ititj = iXiXj.translation();

    Matrix<double, 6, 6> Ji, Jj;

    jacobian_3d_qman ( Ji, Jj,
              izR(0,0), izR(0,1), izR(0,2), izt(0),
              izR(1,0), izR(1,1), izR(1,2), izt(1),
              izR(2,0), izR(2,1), izR(2,2), izt(2),
              iRiRj(0,0), iRiRj(0,1), iRiRj(0,2), ititj(0),
              iRiRj(1,0), iRiRj(1,1), iRiRj(1,2), ititj(1),
              iRiRj(2,0), iRiRj(2,1), iRiRj(2,2), ititj(2));

    _jacobianOplus[0] = Ji;
    _jacobianOplus[1] = Jj;


    _jacobianOplus[0]*=vSwitch->estimate();
    _jacobianOplus[1]*=vSwitch->estimate();


    // derivative w.r.t switch vertex
    _jacobianOplus[2].setZero();

    g2o::SE3Quat delta = measurement().inverse() * iXiXj;
    ErrorVector error;
    error.head<3>() = delta.translation();
    // The analytic Jacobians assume the error in this special form (w beeing positive)
    if (delta.rotation().w() < 0.)
      error.tail<3>() =  - delta.rotation().vec();
    else
      error.tail<3>() =  delta.rotation().vec();

    _jacobianOplus[2] = error * vSwitch->gradient();

}


// ================================================
void EdgeSE3Switchable::computeError()
{
    const VertexSwitchLinear* v3 = static_cast<const VertexSwitchLinear*>(_vertices[2]);
    const g2o::VertexSE3* v1 = dynamic_cast<const g2o::VertexSE3*>(_vertices[0]);
    const g2o::VertexSE3* v2 = dynamic_cast<const g2o::VertexSE3*>(_vertices[1]);

    Eigen::Isometry3d iXiXj_iso = v1->estimate().inverse() * v2->estimate();
    g2o::SE3Quat iXiXj = g2o::SE3Quat(iXiXj_iso.rotation(), iXiXj_iso.translation());
    g2o::SE3Quat delta = measurement().inverse() * iXiXj;
    _error.head<3>() = delta.translation()* v3->estimate();
    
    // The analytic Jacobians assume the error in this special form (w beeing positive)
    if (delta.rotation().w() < 0.)
      _error.tail<3>() =  - delta.rotation().vec()* v3->estimate();
    else
      _error.tail<3>() =  delta.rotation().vec()* v3->estimate();

    
    return;
}
