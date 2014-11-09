#ifndef DAVIDSON_NEW_H
#define DAVIDSON_NEW_H
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "common.h"
#include "mointegrals.h"
#include "cellinfo.h"
#include "ao_ints.h"
void Davidson_Gamma_MPI(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoInts,
    VMatrixXd& inevecs,
    VMatrixXd& inevals,
    const int numberOfEvals,
    const int itol
);

void MP2_Gamma_MPI(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoInts,
    VMatrixXd& inevecs,
    VMatrixXd& inevals,
    const int itol
);

#endif
