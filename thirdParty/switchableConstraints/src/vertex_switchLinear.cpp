#include "vertex_switchLinear.hpp"
#include <iostream>

using namespace std;

VertexSwitchLinear::VertexSwitchLinear() : BaseVertex<1, double>()
{
  _x=0;
  _estimate=_x;
}

bool VertexSwitchLinear::read(std::istream& is)
{
  return true;
}

bool VertexSwitchLinear::write(std::ostream& os) const
{
  return os.good();
}

void VertexSwitchLinear::setToOriginImpl()
{
  _x=0;
  _estimate=_x;

  return;
}


bool VertexSwitchLinear::setEstimateDataImpl(const number_t* et)
{
  _x=(*et);
  _estimate=_x;

  return true;
}


void VertexSwitchLinear::oplusImpl(const number_t* update)
{
  _x += update[0];

  if (_x<0) _x=0;
  if (_x>1) _x=1;

  _estimate=_x;

  return;
}
