#include "edge_se2Switchable.hpp"
#include "vertex_switchLinear.hpp"

using namespace std;
using namespace Eigen;

// ================================================
EdgeSE2Switchable::EdgeSE2Switchable() : g2o::BaseMultiEdge<3, g2o::SE2>()
{
  resize(3);
  _jacobianOplus[0].resize(3,3);    
  _jacobianOplus[1].resize(3,3);
  _jacobianOplus[2].resize(3,1); 

}
// ================================================
bool EdgeSE2Switchable::read(std::istream& is)
  {
    return true;
  }
// ================================================
bool EdgeSE2Switchable::write(std::ostream& os) const
{
    return os.good();
}

// ================================================
void EdgeSE2Switchable::linearizeOplus()
{
   
    const g2o::VertexSE2* vi = static_cast<const g2o::VertexSE2*>(_vertices[0]);
    const g2o::VertexSE2* vj = static_cast<const g2o::VertexSE2*>(_vertices[1]);
    const VertexSwitchLinear* vSwitch = static_cast<const VertexSwitchLinear*>(_vertices[2]);

    double thetai = vi->estimate().rotation().angle();

    Vector2d dt = vj->estimate().translation() - vi->estimate().translation();
    double si=sin(thetai), ci=cos(thetai);

    _jacobianOplus[0](0, 0) = -ci; _jacobianOplus[0](0, 1) = -si; _jacobianOplus[0](0, 2) = -si*dt.x()+ci*dt.y();
    _jacobianOplus[0](1, 0) =  si; _jacobianOplus[0](1, 1) = -ci; _jacobianOplus[0](1, 2) = -ci*dt.x()-si*dt.y();
    _jacobianOplus[0](2, 0) =  0;  _jacobianOplus[0](2, 1) = 0;   _jacobianOplus[0](2, 2) = -1;

    _jacobianOplus[1](0, 0) = ci; _jacobianOplus[1](0, 1)= si; _jacobianOplus[1](0, 2)= 0;
    _jacobianOplus[1](1, 0) =-si; _jacobianOplus[1](1, 1)= ci; _jacobianOplus[1](1, 2)= 0;
    _jacobianOplus[1](2, 0) = 0;  _jacobianOplus[1](2, 1)= 0;  _jacobianOplus[1](2, 2)= 1;

    const g2o::SE2& rmean = measurement().inverse();
    Matrix3d z = Matrix3d::Zero();
    z.block<2, 2>(0, 0) = rmean.rotation().toRotationMatrix();
    z(2, 2) = 1.;
    _jacobianOplus[0] = z * _jacobianOplus[0];
    _jacobianOplus[1] = z * _jacobianOplus[1];



    _jacobianOplus[0]*=vSwitch->estimate();
    _jacobianOplus[1]*=vSwitch->estimate();


    // derivative w.r.t switch vertex
    _jacobianOplus[2].setZero();
    g2o::SE2 delta = measurement().inverse() * (vi->estimate().inverse()*vj->estimate());
    _jacobianOplus[2] = delta.toVector() * vSwitch->gradient();
  
    return;
}


// ================================================
void EdgeSE2Switchable::computeError()
{
    const g2o::VertexSE2* v1 = static_cast<const g2o::VertexSE2*>(_vertices[0]);
    const g2o::VertexSE2* v2 = static_cast<const g2o::VertexSE2*>(_vertices[1]);
    const VertexSwitchLinear* v3 = static_cast<const VertexSwitchLinear*>(_vertices[2]);
  
    g2o::SE2 delta = measurement().inverse() * (v1->estimate().inverse()*v2->estimate());
    _error = delta.toVector() * v3->estimate();
    
    return;
}