#ifndef GENERATE_TVECS
#define GENERATE_TVECS

#include <vector>
#include "Eigen/Dense" 

void generate_tvecs( 
  int& nvals,
  const Eigen::Matrix3d& Tmat,
  int in_xold, int in_xnew,
  std::vector< double >& xyzvals,
  std::vector< double >& xyz_sq_vals
);

void generate_tvecs(
  int& nvals,
  const std::vector< double >& tmatrix,
  int in_xold, int in_xnew,
  int in_yold, int in_ynew,
  int in_zold, int in_znew,
  double in_xoffset,
  double in_yoffset,
  double in_zoffset,
  std::vector< double >& xyzvals,
  std::vector< double >& xyz_sq_vals
);


#endif
