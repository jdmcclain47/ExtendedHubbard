#ifndef FIND_BOUNDING_SPHERE_H
#define FIND_BOUNDING_SPHERE_H

#include "Eigen/Dense"
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
);

#endif
