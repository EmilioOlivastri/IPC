#include "edge_switchPrior.hpp"

using namespace std;

bool EdgeSwitchPrior::read(std::istream &is)
{
  return true;
}

bool EdgeSwitchPrior::write(std::ostream &os) const
{
  return true;
}

void EdgeSwitchPrior::linearizeOplus()
{
  _jacobianOplusXi[0]=-1.0;
}

void EdgeSwitchPrior::computeError()
{
  const VertexSwitchLinear* s = static_cast<const VertexSwitchLinear*>(_vertices[0]);

  _error[0] = measurement() - s->x();

}

