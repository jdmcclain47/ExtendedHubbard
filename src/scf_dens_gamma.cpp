#include <vector>
#include <iostream>
#include "scf.h"
#include "Eigen/Dense"
#include "cellinfo.h"
#include "common.h"

using namespace std;
using namespace Eigen;

void SCF::build_dens_gamma( UnitCell& UCell, SuperCell& SCell ){
   MatrixXd matrval;
   matrval.resize(nmo_scell,nmo_scell);
   matrval*=0.0;

   vector<double> eigenval_list;
   eigenval_list.resize( nmo_scell );
   for(int j=0;j<nmo_scell;j++){
      eigenval_list[ j ] =  e_vals.irrep(0)(j,0) ;
   }
   sort(eigenval_list.begin(),eigenval_list.end());

   /* With some systems such as tpa, we may have a problem in that a restricted HF approach
    * is not suitable, as the HOMO level may have unfilled orbitals.  One can see this is the case
    * for a system such as a 2n cyclic carbon system, where each of the HOMO's are occupied by 1
    * electron.  The following is a quick check..
    */

   double evaltol = 1e-10;
   if(abs(eigenval_list[nocc] - eigenval_list[nocc - 1]) < evaltol){
      cout << "***********************************************" << endl;
      cout << " HOMO has near degeneracy and RHF may not be " << endl;
      cout << "suitable. Increasing number of cells considered " << endl;
      cout << "by 1 may remove this error..." << endl;
      cout << "***********************************************" << endl;
   }

   vector< int > occupancy;
   occupancy.resize( nmo_scell );
   int ecounter = 0;
   for( int j = 0; j < nmo_scell; j++ ){
      if( (e_vals.irrep(0)(j,0) < eigenval_list[ nocc ]) && (ecounter < (2*nocc)) ){
         occupancy[ j ] = 2;
         ecounter += 2;
      }else{
         occupancy[ j ] = 0;
      }
   }

   double sum;
   for( int i = 0; i < nmo_scell; i++ ){
      for( int j = 0; j <= i ; j++ ){
         sum = 0.0;
         for( int k = 0; k < nmo_scell; k++ ){
            sum += occupancy[ k ] * ( e_vecsXd.irrep( 0 )( i, k ) * e_vecsXd.irrep( 0 )( j, k ) );
         }
         SCell_DensXd.irrep( 0 )( i, j ) = sum;
         SCell_DensXd.irrep( 0 )( j, i ) = sum;
      }
   }


   for(int i=0;i<SCell.total_number_of_cells;i++){
      DensXd.irrep(i) = SCell_DensXd.irrep(0).block(0,i*nmo_ucell,nmo_ucell,nmo_ucell);
   }
}


