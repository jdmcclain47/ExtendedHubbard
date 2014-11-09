#ifndef MO_INTEGRALS_H
#define MO_INTEGRALS_H

#include <vector>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "cellinfo.h"
#include "common.h"
#include "ao_ints.h"
#include "Eigen/Dense"
#include <time.h>

class moIntegralFactory {
    VMatrixXd evals;
    VMatrixXd evecsXd;
    VMatrixXcd evecsXcd;
    int nkpt, ntrans, nmo_ucell, nmo_scell;
    Eigen::MatrixXd kernel_matr;
    Eigen::MatrixXd one_body_matr;
    Eigen::MatrixXd coulomb_matr;
    Eigen::VectorXd gamma_contr;
    Eigen::VectorXd gamma_ij;
    Eigen::VectorXd gamma_kl;

    Eigen::VectorXd gamma_i;
    Eigen::VectorXd gamma_j;
    Eigen::VectorXd gamma_k;
    Eigen::VectorXd gamma_l;

    double gamma_transform( const int& imo, const int& jmo, const int& kmo, const int& lmo, const bool& create_new_kl ); 
    double quarter_gamma_transform( 
      const int& imo, const int& jmo, const int& kmo, const int& lmo, 
      const bool& create_new_i, const bool& create_new_j, const bool& create_new_k, const bool& create_new_l
    ); 

    public :
      int get_nmo_scell(){ return nmo_scell; };
      void Init(
          UnitCell& UCell,
          SuperCell& SCell,
          aoIntegralFactory& aoints,
          const PBC_CLASS& pbc,
          VMatrixXd& inevals,
          VMatrixXd& inevecsXd,
          VMatrixXcd& inevecsXcd
      );
      Eigen::MatrixXd get_kernel_matr(){ return kernel_matr; };
      Eigen::MatrixXd get_evecsXd(){ return evecsXd.irrep( 0 ); };
      std::vector< double > full_gamma_transform();
      std::vector< double > full_gamma_transform(
          int pst, int pend, int qst, int qend, int rst, int rend, int sst, int send
      );
      void check_gamma_exchange_energy( std::vector< double >& gamma_ints, const double& ref_energy );
      void gamma_moint_to_vector( std::vector< double >& dble_arr, std::vector< size_t >& indx_arr, size_t& arr_size, 
                                  int pstart, int pend, int qstart, int qend, int rstart, int rend, int tstart, int tend, int itol );
      void write_gamma_moint_to_file( const char* ofilename );
      void write_gamma_moint_to_fcidump( const char* ofilename, int tol );
      void write_gamma_mointb_ind_p( 
          std::string ofilename,
          int tol
      );
      std::vector< double > read_gamma_moints_from_file( const char* ifilename );
};

double quarter_transform__gamma_point(
    const Eigen::VectorXd& gamma_i, 
    const Eigen::VectorXd& gamma_j, 
    const Eigen::VectorXd& gamma_k, 
    const Eigen::VectorXd& gamma_l, 
    Eigen::VectorXd& gamma_contr1,
    Eigen::VectorXd& gamma_contr2,
    Eigen::VectorXd& gamma_contr3,
    const Eigen::MatrixXd& kernel_matr, 
    const int nmo_scell,
    const bool& imo_is_new, 
    const bool& jmo_is_new, 
    const bool& kmo_is_new, 
    const bool& lmo_is_new
);

void gamma_moint_to_vector( 
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    const Eigen::MatrixXd& kernel_matr,
    const Eigen::MatrixXd& evecsXd,
    const int nmo_scell,
    size_t& arr_size,
    int pstart, int pend,
    int qstart, int qend,
    int rstart, int rend,
    int sstart, int send,
    int itol
);

#endif;
