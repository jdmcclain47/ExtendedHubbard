#ifndef DAVIDSON_H
#define DAVIDSON_H
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "common.h"
#include "mointegrals.h"
#include "cellinfo.h"
struct eigenpairs {
   Eigen::VectorXd evals;
   Eigen::MatrixXd evecs;
};

struct eigenpairsXcd {
   Eigen::VectorXd  evals;
   Eigen::MatrixXcd evecs;
};

class Davidson {

private:
  std::vector< double > dble_arr; 
  std::vector< size_t > indx_arr; 
  size_t size_arr;
  std::vector< double > dble_arr2; 
  std::vector< size_t > indx_arr2; 
  size_t size_arr2;
  size_t opt_size_arr;
  bool optimize_size_arr;

   Eigen::MatrixXcd c_hmatr;
   Eigen::MatrixXd r_hmatr;
   VMatrixXd evals;
   Eigen::MatrixXcd c_evecs;
   Eigen::MatrixXd r_evecs;  
   Eigen::VectorXd Hdiag;
   Eigen::VectorXcd HdiagXcd;
   int ckpoint;

   std::vector< double > moints;
   std::vector< std::complex< double > > kmoints;
   UnitCell UCell;
   SuperCell SCell;

   bool full_moint_set;
   double tol;
   int max_sub;
   int nocc, norb;
   int cidim;
   int _nNeighbors;
   moIntegralFactory cMOints;

public:
   void Init(
      UnitCell&,
      SuperCell&,
      moIntegralFactory&
   );
   void set_evals_and_moints(
      VMatrixXd& evals,
      std::vector< double >& moints
   );
   void set_evals( VMatrixXd& evals );
   void set_evals_and_moints(
      VMatrixXd& evals,
      std::vector< std::complex< double > >& moints
   );
   void DavidsonCIS(int sub_size,int numberOfEvals, int tol, std::string mointfile );
   void DavidsonCIS_kpoint( int sub_size, int numberOfEvals, int tol, bool write_eigenvalues, const char* outfile);
   void DavidsonCIS_kpoint( int neval, int tol ){ DavidsonCIS_kpoint( 30, neval, tol, false, ""); };
   void DavidsonCIS(int sub_size,int numberOfEvals, int nNeighbors, int tol);
   Eigen::VectorXcd fill_hrvecXcd( Eigen::VectorXcd& cvec );
   Eigen::VectorXd Fill_Hrvec(Eigen::VectorXd _cvec, std::string mointfile);
   void CheckerCIS();
   void CheckerkCIS(int inkpoint, bool write_eigenvalues, std::ofstream& outstream );
};

#endif

