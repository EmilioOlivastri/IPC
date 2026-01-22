#pragma once

#include "g2o/core/g2o_core_api.h"
#include "g2o/core/robust_kernel.h"

namespace g2o {

class G2O_CORE_API RobustKernelIPC : public RobustKernel {
 public:
  void robustify(number_t error, Vector3& rho) const;
  void setScale(number_t scale);

 protected:
  number_t _scale;
};

class G2O_CORE_API RobustKernelIPCGatedConstant : public RobustKernel {
 public:
  void robustify(number_t error, Vector3& rho) const;
  void setScale(number_t scale);

 protected:
  number_t _scale;
};

class G2O_CORE_API RobustKernelIPCQuad : public RobustKernel {
 public:
  virtual void robustify(number_t e2, Vector3& rho) const;
};

class G2O_CORE_API RobustKernelIPCGatedQuad : public RobustKernel {
 public:
  virtual void robustify(number_t e2, Vector3& rho) const;
};

}  // end namespace g2o
