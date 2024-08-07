#include <Eigen/Core>
#include <Eigen/Dense>

void jacobian_3d_qman (Eigen::Matrix<double, 6, 6> &  Ji , Eigen::Matrix<double, 6, 6> &  Jj,
    const double&  z11 , const double&  z12 , const double&  z13 , const double&  z14 ,
  const double&  z21 , const double&  z22 , const double&  z23 , const double&  z24 ,
  const double&  z31 , const double&  z32 , const double&  z33 , const double&  z34 ,
  const double&  xab11 , const double&  xab12 , const double&  xab13 , const double&  xab14 ,
  const double&  xab21 , const double&  xab22 , const double&  xab23 , const double&  xab24 ,
  const double&  xab31 , const double&  xab32 , const double&  xab33 , const double&  xab34 );