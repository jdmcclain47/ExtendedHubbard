#ifndef MP2_H
#define MP2_H

#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include "Eigen/Dense"
#include "common.h"
#include "cellinfo.h"
#include "mointegrals.h"

double mp2_gamma_full_moints(
    UnitCell& UCell,
    SuperCell& SCell,
    VMatrixXd& evals,
    std::vector< double >& moints
);

double mp2_gamma_ind_p(
  UnitCell& UCell,
  SuperCell& SCell,
  VMatrixXd& evals,
  std::string mointfile,
  bool read_in_from_file,
  moIntegralFactory& moint_class,
  int argc, char *argv[]
);

#endif
