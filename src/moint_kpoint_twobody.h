#ifndef TWO_BODY_K
#define TWO_BODY_K

#include <vector>
#include <complex>
#include "cellinfo.h"
#include "common.h"
#include "ao_ints.h"
#include "Eigen/Dense"

class moint_driver_k{

    double self_interaction_energy;
    int nkpt;
    int ntrans;
    int kernel;
    int nmo_ucell, nmo_scell;
    VMatrixXcd evecsXcd;
    VMatrixXd evals;
    Eigen::MatrixXd kernel_matr; 

    public:
      void kpoint_init( 
             UnitCell& inUCell,
             SuperCell& inSCell,
             aoIntegralFactory& aoints,
             VMatrixXd& inevals,
             VMatrixXcd& inevecsXcd
           );
      std::vector< std::complex< double > > create_kpoint_2e_ints(
             UnitCell& inUCell,
             SuperCell& inSCell,
             aoIntegralFactory& aoints,
             VMatrixXcd& inevecsXcd
      );
      std::vector< std::complex< double > > create_kpoint_2e_ints_fast_FULL(
             UnitCell& inUCell,
             SuperCell& inSCell,
             aoIntegralFactory& aoints,
             VMatrixXcd& inevecsXcd
      );
      std::vector< std::complex< double > > create_kpoint_2e_ints_fast(
          UnitCell& inUCell,
          SuperCell& inSCell,
          aoIntegralFactory& aoints,
          VMatrixXcd& inevecsXcd,
          int pstart,
          int plast,
          int qstart,
          int qlast,
          int rstart,
          int rlast,
          int sstart,
          int slast
      );
      void check_exchange_energy( 
          UnitCell& inUCell,
          SuperCell& inSCell,
          aoIntegralFactory& aoints,
          std::vector< std::complex< double > >& moints 
      );
      void write_out_moints( std::vector< std::complex< double > > kmoints, const char* ofilename );
      std::vector< std::complex< double > > read_in_moints( const char* ifilename );
      double get_kernel( int mu, int nu );
      void set_kernel( int ktype );

};

#endif
