#ifndef EWALD_H
#define EWALD_H

#include "Eigen/Dense"
#include <vector>

Eigen::Matrix3d create_Gmat( Eigen::Matrix3d& Tmat );

double get_alpha( Eigen::Matrix3d& Tmat );

void setalpha( double alpha_in );

double energy_ewald_converge_FAST( 
   Eigen::Matrix3d& Tvec, 
   Eigen::Matrix3d& Gvec, 
   std::vector< Eigen::Vector3d >& coord,
   std::vector< double >& charge,
   int itol
);

double potential_ewald_converge_FAST( 
   Eigen::Vector3d& pos, 
   Eigen::Matrix3d& Tvec, 
   Eigen::Matrix3d& Gvec, 
   std::vector< Eigen::Vector3d >& coord, 
   std::vector< double >& charge,
   int itol
);

double potential_ewald_converge_FAST( 
   Eigen::Vector3d& pos, 
   Eigen::Matrix3d& Tvec, 
   Eigen::Matrix3d& Gvec, 
   Eigen::Vector3d& coord, 
   std::vector< double >& charge,
   int itol
);

double self_potential_ewald_converge_FAST(
   Eigen::Matrix3d& Tvec, 
   Eigen::Matrix3d& Gvec, 
   int itol
);

double energy_ewald( 
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec,
    Eigen::Matrix< double, Eigen::Dynamic, 3 >& coord,
    std::vector<double>& charge
);

double energy_ewald(
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec,
    std::vector< Eigen::Vector3d >& coord,
    std::vector<double>& charge
);

double energy_ewald_converge(
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec,
    std::vector< Eigen::Vector3d >& coord,
    std::vector<double>& charge,
    int itol
);

double potential_ewald( 
    Eigen::Vector3d pos, 
    Eigen::Matrix3d& Tvec, 
    Eigen::Matrix3d& Gvec, 
    Eigen::Matrix< double, Eigen::Dynamic, 3 >& coord, 
    std::vector<double>& charge
);

double potential_ewald( 
    Eigen::Vector3d pos, 
    Eigen::Matrix3d& Tvec, 
    Eigen::Matrix3d& Gvec, 
    std::vector< Eigen::Vector3d >& coord,
    std::vector<double>& charge
);

double potential_ewald( Eigen::VectorXd pos, Eigen::Matrix3d& Tvec, Eigen::Matrix3d& Gvec, Eigen::Matrix< double, Eigen::Dynamic, 3 >& coord, std::vector<double>& charge );
double potential_ewald( Eigen::Vector3d pos, Eigen::Matrix3d& Tvec, Eigen::Matrix3d& Gvec, Eigen::Vector3d& coord, std::vector< double >& charge ); 
double self_potential_ewald( Eigen::Matrix3d& Tvec, Eigen::Matrix3d& Gvec );
double self_potential_ewald_converge( Eigen::Matrix3d& Tvec, Eigen::Matrix3d& Gvec, int itol );
double potential_ewald_converge( Eigen::Vector3d pos, Eigen::Matrix3d& Tvec, Eigen::Matrix3d& Gvec, 
    Eigen::Vector3d& coord, std::vector< double >& charge, int itol );

double potential_ewald__2D( 
    Eigen::VectorXd pos, 
    Eigen::Matrix3d& Tvec, 
    Eigen::Matrix3d& Gvec, 
    Eigen::Matrix< double, Eigen::Dynamic, 3 >& coord, 
    std::vector< double >& charge
);
double self_potential_ewald__2D(
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec
);

#endif
