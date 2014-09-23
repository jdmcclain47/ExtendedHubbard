#ifndef DIIS_H
#define DIIS_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include "Eigen/Dense" 
#include <fstream>
#include <vector>

class diis {
   int nmax;
   int nblocks;
   std::vector< Eigen::MatrixXd > focklist;
   std::vector< Eigen::MatrixXd > errlist;

   public :
      void restart();
      void Init( int nmax_, int nblocks_ );
      void use_diis( Eigen::MatrixXd* outFock, Eigen::MatrixXd inFock, Eigen::MatrixXd inDens );
};

#endif
