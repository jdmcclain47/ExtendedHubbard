#include "Eigen/Dense"
#include "find_bounding_sphere.h"
#include <math.h>

// For a given distance and translation vectors, this finds how many translation vectors
// can fit a sphere of that given distance within its boundaries.
//
// a.k.a. finds how far you need to go out in each translation direction to look for all
// orbitals/atoms within a given distance

void findBoundingSphere(
    const Eigen::Matrix3d& Tmat,
    double rcut,
    int& nx, 
    int& ny, 
    int& nz
){
    double rtlowesteigen;
    Eigen::MatrixXd TxTtrans = Tmat * ( Tmat.transpose() );
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > umin( TxTtrans );
    Eigen::VectorXd eig = umin.eigenvalues();
    rtlowesteigen = eig.minCoeff();
    rtlowesteigen = sqrt( rtlowesteigen );

    nx = (int)std::ceil( rcut / rtlowesteigen );
    ny = (int)std::ceil( rcut / rtlowesteigen );
    nz = (int)std::ceil( rcut / rtlowesteigen );
    return;
}
