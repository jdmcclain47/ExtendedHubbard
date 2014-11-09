#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include "Eigen/Dense"
#include "davidson.h"
#include <time.h>
#include "binary_search.h"
#include "moint_read_file.h"
#include "mointegrals.h"
#include "common.h"
#include "sorter.h"
#include "cellinfo.h"

using namespace std;
using namespace Eigen;

static size_t get_ijkl( int i, int j, int k, int l ){
    size_t pq, rs, index;
    pq = ( i > j ) ? (i*(i+1))/2+j : (j*(j+1))/2+i;
    rs = ( k > l ) ? (k*(k+1))/2+l : (l*(l+1))/2+k;
    index = ( pq > rs ) ? (pq*(pq+1))/2+rs : (rs*(rs+1))/2+pq; 
    return index;
}

void Davidson::DavidsonCIS(
    int sub_size,
    int numberOfEvals_,
    int tol_,
    std::string mointfile
){

  /* Setting everything up for MOINTEGRALS */

  dble_arr.resize( 1000000 );
  indx_arr.resize( 1000000 );
  dble_arr2.resize( 1000000 );
  indx_arr2.resize( 1000000 );
  opt_size_arr = 0;
  optimize_size_arr = true;
//  if( !full_moint_set ){
//     std::cout << "READING FROM MULTIPLE MOINT FILES" << std::endl;
//     for( int i = 0; i < SCell.nao; ++i ){
//       read_gamma_mointb_ind_p( mointfile, i, dble_arr, indx_arr, size_arr, optimize_size_arr );
//       if( size_arr > opt_size_arr ) opt_size_arr = size_arr;
//     }
//     std::cout << "OPTIMAL SIZE OF ARRAYS " << opt_size_arr << std::endl;
//     optimize_size_arr = false;
//     dble_arr.resize( opt_size_arr );
//     indx_arr.resize( opt_size_arr );
//
//     dble_arr2.resize( opt_size_arr );
//     indx_arr2.resize( opt_size_arr );
///*
//     for( int i = 0; i < moints.get_nmo_scell(); ++i ){
//       read_gamma_mointb_ind_p( mointfile, i, dble_arr, indx_arr, size_arr, optimize_size_arr );
//       int INDEX_TO_BE_FOUND = 1673535;
//       int mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
//       if( mo_index_in_arr >= 0 ) std::cout << "MOINTEGRAL " << INDEX_TO_BE_FOUND << " = " << dble_arr[ mo_index_in_arr ] << std::endl;
//     }
//*/
//  }

  nocc = (int)(SCell.nao/2);
  norb = SCell.nao;
  cidim = nocc * nocc;

   tol = pow( 10., -1. * tol_ );

   int psia,psib,psir,psis;
   int pq, rs, index1, index2;
   double cival;
   int numberOfEvals = numberOfEvals_;
   int index_ci;

   //For more than one eigenvalue, it more just checks that the other eigenvalues are
   //   eigenvalues, rather than "lowest" eigenvalues.  Might need to change which
   //   current_eval you are trying to get close to

   if( numberOfEvals > cidim ){
      cout << "***********************************************************************" << endl;
      cout << "Trying to find more eigenvalues than what's possible in this subspace.." << endl;
      cout << "  Will instead find max number of eigenvalues.                         " << endl;
      cout << "***********************************************************************" << endl;
      numberOfEvals = cidim;
   }

   // sub_size to control max subspace size is currently not implemented.. just stops if subspace
   // has reached sub_size dimension (rather than replacing some of the old vectors with new ones)

   // Making the diagonal Matrix
   Hdiag.resize(cidim);
   double moint1, moint2;
   int mo_index_in_arr;
   size_t INDEX_TO_BE_FOUND;
   std::cout << "CREATING DIAGONAL MATRIX..." << std::endl;
   clock_t t = clock();
   int old_psir = -1;
   size_t lower_bound1, lower_bound2;
   lower_bound1 = 0;
   lower_bound2 = 0;

//void moIntegralFactory::vec_gamma_mointb_ind_p( 
//    const int& p1,
//    const int& p2,
//    std::vector< double >& dble_arr,
//    std::vector< size_t >& indx_arr,
//    size_t& size_arr,
//    int tol
//){

   for( int psir = nocc; psir < norb; psir++ ){
   for( int psia = 0; psia < nocc; psia++ ){
//
//   for( int i = 0; i < cidim; i++ ){
//      psia = i - (int)(i / nocc) * nocc;
//      psir = ( i / nocc ) + nocc;
      psib = psia;
      psis = psir;

      // ... psir > psia
      pq = ( psir * ( psir + 1 ) ) / 2 + psia;
      rs = ( psis * ( psis + 1 ) ) / 2 + psib;
      index1 = ( pq * ( pq + 1 ) ) / 2 + rs; 

      pq = ( psir * ( psir + 1 ) ) / 2 + psis;
      rs = ( psia * ( psia + 1 ) ) / 2 + psib;
      index2 = ( pq * ( pq + 1 ) ) / 2 + rs; 

      // ... diagonal for RHF is just 
      //        eps(r) - eps(a) + 2*(ra|bs) - (rs|ba)
      cival  = evals.irrep( 0 )( psir, 0 ) - evals.irrep( 0 )( psia, 0 );

      if( full_moint_set ){
        cival += 2. * moints[ index1 ] - 1. * moints[ index2 ];
      }else{
        if( old_psir != psir ){
          int psirp1 = psir+1;
          cMOints.vec_gamma_mointb_ind_p( psir, psirp1, dble_arr, indx_arr, size_arr, -10 );
        }

        moint1 = 0.0;
        INDEX_TO_BE_FOUND = index1;
        mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
        if( mo_index_in_arr >= 0 ){
          moint1 = dble_arr[ mo_index_in_arr ];
          lower_bound1 = mo_index_in_arr;
        }
        moint2 = 0.0;
        INDEX_TO_BE_FOUND = index2;
        mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
        if( mo_index_in_arr >= 0 ){
          moint2 = dble_arr[ mo_index_in_arr ];
          lower_bound1 = mo_index_in_arr;
        }
        cival += 2. * moint1 - 1. * moint2;
        
        old_psir = psir;
      }

      index_ci = (psir-nocc) * nocc + psia;
      
      Hdiag( index_ci ) = cival;
      Hdiag( index_ci ) = 0.0;
   }
   }
   //cout << "HDIAGONAL " << endl;
   //cout << Hdiag << endl;
   //cout << "--------------------" << endl;

   t = clock() - t;
   std::cout << " - done in " << t / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
   int cdim = cidim;
   int CurrentSize;

   vector<VectorXd> r_evecs;
   VectorXd tempvec = VectorXd::Ones(cdim);
   r_evecs.push_back(tempvec);
   CurrentSize = 1;

   VectorXd residual;
   residual.resize(cdim);
   VectorXd corr_vec;
   vector<VectorXd> Ar_evecs;
   eigenpairs _eigenpairs;

   corr_vec.resize(cdim);

   MatrixXd identitymatr;

   int CorrectedSize=CurrentSize;
   int OldSize=CurrentSize;

   // initialization of first vector and Hr vector
   std::cout << "INITIALIZING FIRST VECTOR..." << std::endl;
   t = clock();
   for(int i=0;i<CurrentSize;i++){
      if(i==0){
         r_evecs[0]=r_evecs[0].normalized();
         Ar_evecs.push_back(Fill_Hrvec(r_evecs[0], mointfile));
         //cout << "R VEC " << endl;
         //cout << r_evecs[ i ].transpose() << endl;
         //cout << "FIRST HR VEC " << endl;
         //cout << Ar_evecs[ i ].transpose() << endl;
      }else{
         for(int j=0;j<i;j++){
            r_evecs[i] -=  (r_evecs[j] * r_evecs[i].transpose() * r_evecs[j] ) ;
         }
         r_evecs[i]=r_evecs[i].normalized();
         Ar_evecs.push_back(Fill_Hrvec(r_evecs[i], mointfile));
         //cout << "R VEC " << endl;
         //cout << r_evecs[ i ].transpose() << endl;
         //cout << "HR VEC " << endl;
         //cout << Ar_evecs[ i ].transpose() << endl;
      }
   }
   t = clock() - t;
   std::cout << " - done in " << t / (double)CLOCKS_PER_SEC << " seconds." << std::endl;

   int convergedevals_=0;
   int converged=0;

   // Setting up the subspace matrix
   MatrixXd OldSubspaceMatr;
   OldSubspaceMatr.resize(OldSize,OldSize);
   for(int i=0;i<OldSize;i++){
      for(int j=0;j<=i;j++){
         OldSubspaceMatr(i,j) = (r_evecs[i].transpose()) * Ar_evecs[j];
         OldSubspaceMatr(j,i) = OldSubspaceMatr(i,j);
      }
   }
   MatrixXd SubspaceMatr;

   //starting the convergence cycles
   while(converged==0){

      //Orthonormalizing current space of vectors
      //orthogonalizing vectors if new subspace is larger than old one
      if(OldSize!=CurrentSize){	
         for(int i=OldSize;i<CurrentSize;i++){
	     for(int j=0;j<i;j++){
	        r_evecs[i] -=  (r_evecs[j] * r_evecs[i].transpose() * r_evecs[j] ) ;
	     }
	     r_evecs[i]=r_evecs[i].normalized();
	     Ar_evecs.push_back(Fill_Hrvec(r_evecs[i], mointfile));
         //cout << "R VEC " << endl;
         //cout << r_evecs[ i ].transpose() << endl;
         //cout << "HR VEC " << endl;
         //cout << Ar_evecs[ i ].transpose() << endl;
         }
      } 

      /* Now we have the r vectors and the Hr vectors.. Now we form the r^T H r subspace matrix
       * noting that for all the r's and Hr's that were kept in the previous iteration, we only
       * need to calculate the r^T H r elements for the newly added r's and Hr's
       */ 
      // This is our current Subspace we will be diagonalizing
      SubspaceMatr.resize(CurrentSize,CurrentSize);
      // re-using parts of the old subspace matrix and then calculating new subspace matrix
      SubspaceMatr.block(0,0,OldSize,OldSize) = OldSubspaceMatr;
      for(int i=0;i<CurrentSize;i++){
         for(int j=OldSize;j<CurrentSize;j++){
            SubspaceMatr(i,j) = (r_evecs[i].transpose()) * Ar_evecs[j];
            SubspaceMatr(j,i) = SubspaceMatr(i,j);
         }
      }
      OldSubspaceMatr.resize(CurrentSize,CurrentSize);
      OldSubspaceMatr = SubspaceMatr;
      OldSize = CurrentSize;

      SelfAdjointEigenSolver<MatrixXd> eigensolver(SubspaceMatr);
      _eigenpairs.evals = eigensolver.eigenvalues();
      _eigenpairs.evecs = eigensolver.eigenvectors();

      /* Now calculating the residual vector, finding the eigenfunction in the original basis
       * and calculating |A*v - e*v|.  Finding lowest nth eigenvalue, so use corresponding
       * eigenvector in the projected subspace.
       */
      residual*=0;
      for(int i=0;i<cdim;i++){
         for(int j=0;j<CurrentSize;j++){
            residual(i)+=Ar_evecs[j](i) * _eigenpairs.evecs(j,convergedevals_);
            residual(i)-=_eigenpairs.evals(convergedevals_) * r_evecs[j](i) * _eigenpairs.evecs(j,convergedevals_);
         }
      }

      //Computing the correcting vector and seeing if we are going to look for a better/more eigenvalues 
      if( (residual.blueNorm() ) > tol){
         //Check to see if we can increase size of subspace
         if(CurrentSize == (sub_size-1)){
            printf("Could not find lowest [%d] eigenvalues with MaxSubspaceSize[%d]\n",numberOfEvals,sub_size);
            converged = 1;
         }
         else{
         //Increase size of subspace, finding correction vector with orthogonal component to current subspace
            for(int i=0;i<cdim;i++){
               corr_vec(i) = -1.0/(Hdiag(i) - _eigenpairs.evals(convergedevals_))*residual(i);
            }
            r_evecs.push_back(corr_vec);
            CorrectedSize++;
            printf("Current eval [%d] : %20.16e   Diff : %20.16e \n",convergedevals_,_eigenpairs.evals(convergedevals_), residual.blueNorm() );
            fflush(stdout);
         }
      }

      //If converged, see whether to go on to a new eigenvalue if it was specified
      if( (residual.blueNorm() ) < tol){
         if(convergedevals_==(numberOfEvals-1)){
            cout << "lowest eigenvalues (H, eV):" << endl;
            for ( int i = 0; i < numberOfEvals; i++ ){
               printf("CIS-CONVERGED EIGENVALUE %3d : %20.16f, %20.16f \n", (i+1), _eigenpairs.evals(i), ( _eigenpairs.evals(i) * 27.211 ) );
            }
            cout << "converged with dimension : " << CurrentSize << endl;
            fflush(stdout);
            converged=1;
         }
         else{
            convergedevals_++;
            // If our current subspace of size N finds the lowest N eigenvalues up to desired accuracy,
            // need to increase the size of our space done so by adding in a random vector
            if(convergedevals_>=CurrentSize){
               tempvec = VectorXd::Random(cdim);
               r_evecs.push_back(tempvec);
               CorrectedSize++;
            }
         }
      }

      CurrentSize=CorrectedSize;   
   }
}



VectorXd Davidson::Fill_Hrvec(VectorXd _cvec, std::string mointfile){

   VectorXd sigma_vec;
   sigma_vec.resize(cidim);

   int psia,psib,psir,psis;
   int pq1, rs1, pq2, rs2, index1, index2, index3;
   
   double cival;
   int mo_index_in_arr;
   size_t INDEX_TO_BE_FOUND;
   double moint1, moint2;
   int pindex1, old_pindex1;
   int pindex2, old_pindex2;

   old_pindex1 = -1; 
   old_pindex2 = -1; 

   /* for reading a specific mo integral from a file, this sets a lower bound on the location of that integral */
   size_t lower_bound1 = 0;
   size_t lower_bound2 = 0;

   int old_psir, old_psis;
   old_psir = -1;
   old_psis = -1;
   int psirp1;

   int ind_ra, ind_rb, ind_sa, ind_sb;
   for(int psir = nocc; psir < norb; psir++ ){
   psirp1 = psir + 1;

   cMOints.gamma_moint_to_vector( 
      dble_arr, indx_arr, 
      size_arr,
      psir, psirp1, 
      0, norb,
      0, nocc,
      0, norb, 
      -10 
   );
   sorter( dble_arr, indx_arr, size_arr );

   for(int psia = 0; psia < nocc; psia++ ){
      ind_ra = (psir-nocc)*nocc + psia;
      cival=0.0;
      //psia = i - (int)( i/nocc ) * nocc;
      //psir = (int)( i / nocc ) + nocc;
      cival=(evals.irrep( 0 )(psir,0)-evals.irrep( 0 )(psia,0))*_cvec( ind_ra );

      pq1 = ( psir * ( psir + 1 ) ) / 2 + psia;

      lower_bound1 = 0;
      lower_bound2 = 0;

      old_psir = psir;

      for(int psis = nocc; psis < norb; psis++ ){

      for(int psib = 0; psib < nocc; psib++ ){
         ind_sb = (psis-nocc)*nocc + psib;
         ind_sa = (psis-nocc)*nocc + psia;
         ind_rb = (psir-nocc)*nocc + psib;

         // ENTRIES ARE  : eps(r) - eps(a) + 2*(ra|bs) - (rs|ba)

         //if( full_moint_set ){
           //
           //
           // Doing coulomb part
           //
           //
           index1 =  get_ijkl( psir, psia, psib, psis );
           index2 =  get_ijkl( psir, psis, psib, psia );
           moint1 = 0.0;
           INDEX_TO_BE_FOUND = index1;
           mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
           if( mo_index_in_arr >= 0 ){
             moint1 = dble_arr[ mo_index_in_arr ];
             lower_bound1 = mo_index_in_arr;
           }
           cival += 2. * moint1 * _cvec( ind_sb );

           index2 =  get_ijkl( psir, psis, psib, psia );
           moint2 = 0.0;
           INDEX_TO_BE_FOUND = index2;
           mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
           if( mo_index_in_arr >= 0 ){
             moint2 = dble_arr[ mo_index_in_arr ];
             lower_bound2 = mo_index_in_arr;
           }
           cival -= 1. * moint2 * _cvec( ind_sb );
           //if( psir == 3 && psia == 0 ){
           //  printf( "%3d %3d : ADDING (%3d,%3d,%3d,%3d) : %20.16f \n", psir, psia, psir, psia, psib, psis, moint1 );
           //  printf( "%3d %3d : SUBTRA (%3d,%3d,%3d,%3d) : %20.16f \n", psir, psia, psir, psis, psib, psia, moint2 );
           //  cout << "CBJ : " << _cvec( ind_sb ) << endl;
           //}
         //}else{
         //}
      }
      }
   sigma_vec( ind_ra ) = cival;   
   }
   }
return sigma_vec;
}




void Davidson::CheckerCIS(){
    printf( "Doing full cis matrix diagonalizaton... \n" );
    nocc = (int)(SCell.nao/2);
    cidim = nocc * nocc;
    int psia,psib,psir,psis;
    int pq1, rs1, pq2, rs2, index1, index2;
    double cival;

    MatrixXd cis_full;
    cis_full = MatrixXd::Zero( cidim, cidim );
    int orbi, orba, orbj, orbb;
    int ai, jb, ab, ji;
    int aijb, abji;
    //
    // 
    // printing out eigenvalues
    //  
    //  
    //for( int i = 0; i < nocc*2; ++i ) printf( "%20.16f \n", evals.irrep( 0 )( i ) );
    for( int i = 0; i < cidim; ++i ){
        orbi = floor( i / nocc );
        orba = i - floor( int( i / nocc ) ) * nocc + nocc;
        for( int j = 0; j < cidim; ++j ){
            orbj = floor( j / nocc );
            orbb = j - floor( int( j / nocc ) ) * nocc + nocc;
            if( orba == orbb && orbi == orbj )
                cis_full( i, j ) += evals.irrep( 0 )( orba ) - evals.irrep( 0 )( orbi ); 
            // ... we have that    orba > orbi
            ai = ( orba * ( orba + 1 ) ) / 2 + orbi;
            // ... likewise        orbb > orbj
            jb = ( orbb * ( orbb + 1 ) ) / 2 + orbj;
            if( orba > orbb ) ab = ( orba * ( orba + 1 ) ) / 2 + orbb;
            else              ab = ( orbb * ( orbb + 1 ) ) / 2 + orba;

            if( orbi > orbj ) ji = ( orbi * ( orbi + 1 ) ) / 2 + orbj;
            else              ji = ( orbj * ( orbj + 1 ) ) / 2 + orbi;

            if( ai > jb ) aijb = (ai*(ai+1))/2+jb;
            else          aijb = (jb*(jb+1))/2+ai;

            if( ab > ji ) abji = (ab*(ab+1))/2+ji;
            else          abji = (ji*(ji+1))/2+ab;

            printf( "Making < i->a | H - E_0 | j->b > \n" );
            printf( "ijab = %3d, %3d, %3d, %3d \n", orbi, orbj, orba, orbb );
            printf( "eigenvalue contribution : %20.16f \n", cis_full( i, j ) );
            cis_full( i, j ) += 2*moints[ aijb ] - 1.*moints[ abji ];
            printf( "(%3d %3d | %3d %3d ) : %20.16f \n", orba, orbi, orbj, orbb, moints[ aijb ] );
            printf( "(%3d %3d | %3d %3d ) : %20.16f \n", orba, orbb, orbj, orbi, moints[ abji ] );
        }
    }
    //
    // 
    // printing out cis matrix 
    //  
    //  
    // cout << cis_full << endl;
    SelfAdjointEigenSolver< MatrixXd > solver( cis_full );
    cout << "REAL SPACE EIGENVALUES : " << endl;
    for( int i = 0; i < cidim; ++i ) printf( "%20.16f \n", solver.eigenvalues()( i ) );
}
