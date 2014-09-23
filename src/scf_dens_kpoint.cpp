#include <vector>
#include <complex>
#include <iostream>
#include "scf.h"
#include "common.h"
#include "Eigen/Dense"
#include "cellinfo.h"
using namespace std;
using namespace Eigen;

void SCF::build_dens_k( UnitCell& UCell, SuperCell& SCell ){
   MatrixXcd matrval;
   matrval.resize(nmo_ucell,nmo_ucell);
   matrval*=0.0;

   int which_irrep;
   vector<double> eigenval_list;
   for(int i=0;i<nirreps;i++){
      for(int j=0;j<nmo_ucell;j++){
         eigenval_list.push_back( e_vals.irrep(i)(j,0) );
      }
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

   vector< vector< int > > occupancy;
   occupancy.resize( nirreps );
   for( int i = 0; i < nirreps; ++i ){
      occupancy[ i ].resize( nmo_ucell );
   }

   int ecounter = 0;
   for( int i = 0; i < nirreps; ++i ){
      for( int j = 0; j < nmo_ucell; j++ ){
          if( ( e_vals.irrep( i )( j, 0 ) < eigenval_list[ nocc ] ) && (ecounter < (2*nocc)) ){
              occupancy[ i ][ j ] = 2;
              ecounter += 2;
          }else{
              occupancy[ i ][ j ] = 0;
          }
      } 
   }

   complex< double > sum;
   for( int i = 0; i < nirreps; ++i ){
       for( int j = 0; j < nmo_ucell; ++j ){
           // ... note that since we're dealing with the density in reciproal space, we have that it is hermitian
           //     as opposed to the various "blocks" in coordinate space for the density matrix
           for( int k = 0; k <= j; ++k ){
               sum *= 0.0;
               // ... for this irrep, construct the density by summing over all occupied orbitals
               for( int l = 0; l < nmo_ucell; ++l ){
                     real(sum) += occupancy[ i ][ l ] * \
                                  (( real( e_vecsXcd.irrep( i )( j, l )) * real( e_vecsXcd.irrep( i )( k, l ))) + \
                                   ( imag( e_vecsXcd.irrep( i )( j, l )) * imag( e_vecsXcd.irrep( i )( k, l ))));
                     imag(sum) += occupancy[ i ][ l ] * \
                                  (( imag( e_vecsXcd.irrep( i )( j, l )) * real( e_vecsXcd.irrep( i )( k, l ))) - \
                                   ( real( e_vecsXcd.irrep( i )( j, l )) * imag( e_vecsXcd.irrep( i )( k, l ))));
               }
               real( matrval( j, k ) ) = real( sum );
               imag( matrval( j, k ) ) = imag( sum );
               if( k < j ){
                   real( matrval( k, j ) ) = real( sum );
                   imag( matrval( k, j ) ) = - 1. * imag( sum );
               }
           }
       }
       DensXcd.irrep( i ) = matrval;
   }

}
