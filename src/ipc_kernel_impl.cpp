#include "ipc/ipc_kernel_impl.hpp"
//#include "robust_kernel_factory.h"

namespace g2o {

/*----------------------------------------------------------------*/

// ORIGINAL IPC IMPLEMENTATION AS A ROBUST KERNEL
void RobustKernelIPC::setScale(number_t scale) { _scale = scale; }

void RobustKernelIPC::robustify(number_t e2, Vector3& rho) const 
{
  rho[0] = _scale * e2;
  rho[1] = _scale;
  rho[2] = 0.;
}

/*----------------------------------------------------------------*/

// IPC ROBUST KERNEL with Gated Activation
void RobustKernelIPCGatedConstant::setScale(number_t scale) { _scale = scale; }

// The 0.5 should be already included in the error computation
void RobustKernelIPCGatedConstant::robustify(number_t e2, Vector3& rho) const 
{
  number_t dsqr = _delta * _delta;
  if (e2 <= dsqr) // normal quadratic error 
  {  
    rho[0] = e2;
    rho[1] = 1.;
    rho[2] = 0.;
  } 
  else 
  {  
    rho[0] = (_scale * e2 + dsqr * (1. - _scale));
    rho[1] = _scale;
    rho[2] = 0.;
  }
}

/*----------------------------------------------------------------*/

// ALTERNATIVE IPC IMPLEMENTATION AS A ROBUST KERNEL, NO MANUALLY SET SCALING NEEDED
void RobustKernelIPCQuad::robustify(number_t e2, Vector3& rho) const 
{
  rho[0] = 0.5 * e2 * e2;
  rho[1] = e2;
  rho[2] = 1.; // Second derivative not used so not computed
}

/*----------------------------------------------------------------*/

void RobustKernelIPCGatedQuad::robustify(number_t e2, Vector3& rho) const 
{
  number_t dsqr = _delta * _delta;
  if (e2 <= dsqr) // normal quadratic error
  {  
    rho[0] = e2;
    rho[1] = 1.;
    rho[2] = 0.;
  } 
  else // outlier
  {
    rho[0] = 0.5 * e2 * e2 + dsqr * (1. - 0.5 * dsqr);
    rho[1] = e2;
    rho[2] = 1.; // Second derivative not used so not computed
  }
}

/*----------------------------------------------------------------*/

}  // end namespace g2o
