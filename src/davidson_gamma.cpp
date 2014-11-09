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
#include "common.h"
#include "cellinfo.h"

using namespace std;
using namespace Eigen;


void Davidson::DavidsonCIS(
    int sub_size,
    int numberOfEvals_,
    int tol_,
    std::string mointfile
){

  /* Setting everything up for MOINTEGRALS */

  opt_size_arr = 0;
  optimize_size_arr = true;
  if( !full_moint_set ){
     std::cout << "READING FROM MULTIPLE MOINT FILES" << std::endl;
     for( int i = 0; i < SCell.nao; ++i ){
       read_gamma_mointb_ind_p( mointfile, i, dble_arr, indx_arr, size_arr, optimize_size_arr );
       if( size_arr > opt_size_arr ) opt_size_arr = size_arr;
     }
     std::cout << "OPTIMAL SIZE OF ARRAYS " << opt_size_arr << std::endl;
     optimize_size_arr = false;
     dble_arr.resize( opt_size_arr );
     indx_arr.resize( opt_size_arr );

     dble_arr2.resize( opt_size_arr );
     indx_arr2.resize( opt_size_arr );
/*
     for( int i = 0; i < moints.get_nmo_scell(); ++i ){
       read_gamma_mointb_ind_p( mointfile, i, dble_arr, indx_arr, size_arr, optimize_size_arr );
       int INDEX_TO_BE_FOUND = 1673535;
       int mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
       if( mo_index_in_arr >= 0 ) std::cout << "MOINTEGRAL " << INDEX_TO_BE_FOUND << " = " << dble_arr[ mo_index_in_arr ] << std::endl;
     }
*/
  }



   tol = pow( 10., -1. * tol_ );

   int psia,psib,psir,psis;
   int pq, rs, index1, index2;
   double cival;
   int numberOfEvals = numberOfEvals_;
   nocc = (int)(SCell.nao/2);
   cidim = nocc * nocc;

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


   for( int i = 0; i < cidim; i++ ){
      psia = i - (int)(i / nocc) * nocc;
      psir = ( i / nocc ) + nocc;
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
          read_gamma_mointb_ind_p( mointfile, psir, dble_arr, indx_arr, size_arr, optimize_size_arr );
          lower_bound1 = 0;
          lower_bound2 = 0;
        }

        moint1 = 0.0;
        INDEX_TO_BE_FOUND = index1;
        mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, lower_bound1, size_arr-1 );
        if( mo_index_in_arr >= 0 ){
          moint1 = dble_arr[ mo_index_in_arr ];
          lower_bound1 = mo_index_in_arr;
        }


        moint2 = 0.0;
        INDEX_TO_BE_FOUND = index2;
        mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, lower_bound2, size_arr-1 );
        if( mo_index_in_arr >= 0 ){
          moint2 = dble_arr[ mo_index_in_arr ];
          lower_bound2 = mo_index_in_arr;
        }

        cival += 2. * moint1 - 1. * moint2;
        
        old_psir = psir;
      }
      
      Hdiag(i) = cival;
      Hdiag(i) = 0.0;
   }
   t = clock() - t;
   std::cout << " - done in " << t / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
   int cdim = cidim;
   int CurrentSize;


   vector<VectorXd> r_evecs;
   VectorXd tempvec = VectorXd::Random(cdim);
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
      }else{
         for(int j=0;j<i;j++){
            r_evecs[i] -=  (r_evecs[j] * r_evecs[i].transpose() * r_evecs[j] ) ;
         }
         r_evecs[i]=r_evecs[i].normalized();
         Ar_evecs.push_back(Fill_Hrvec(r_evecs[i], mointfile));
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
   int pq1, rs1, pq2, rs2, index1, index2;
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

   for(int i=0;i<cidim;i++){
      cival=0.0;
      psia = i - (int)( i/nocc ) * nocc;
      psir = (int)( i / nocc ) + nocc;
      cival=(evals.irrep( 0 )(psir,0)-evals.irrep( 0 )(psia,0))*_cvec(i);

      pq1 = ( psir * ( psir + 1 ) ) / 2 + psia;

      lower_bound1 = 0;
      lower_bound2 = 0;

      if( !full_moint_set ){
        if( psir != old_psir )
          read_gamma_mointb_ind_p( mointfile, psir, dble_arr, indx_arr, size_arr, optimize_size_arr );
        old_psir = psir;
      }

      for(int j=0;j<cidim;j++){
         psib= j - (int)( j / nocc ) * nocc;
         psis= (int)( j / nocc ) + nocc;

         rs1 = ( psis * ( psis + 1 ) ) / 2 + psib;

         // ENTRIES ARE  : eps(r) - eps(a) + 2*(ra|bs) - (rs|ba)
         //        or... : eps(r) - eps(a) + 2*(ra|sb) - (rs|ba)

         if( pq1 > rs1 ){
             index1 = ( pq1 * ( pq1 + 1 ) ) / 2 + rs1;
             pindex1 = psir;
         }else{
             index1 = ( rs1 * ( rs1 + 1 ) ) / 2 + pq1;
             pindex1 = psis;
         }
         if( psir > psis ){
             pq2 = (psir * ( psir + 1 ) ) / 2 + psis;
             pindex2 = psir;
         }else{
             pq2 = (psis * ( psis + 1 ) ) / 2 + psir;
             pindex2 = psis;
         }
         if( psib > psia ){
             rs2 = (psib * ( psib + 1 ) ) / 2 + psia;
         }else{
             rs2 = (psia * ( psia + 1 ) ) / 2 + psib;
         }
         if( pq2 > rs2 ){
             index2 = ( pq2 * ( pq2 + 1 ) ) / 2 + rs2;
         }else{
             index2 = ( rs2 * ( rs2 + 1 ) ) / 2 + pq2;
         }
         if( full_moint_set ){
           cival += ( 2. * moints[ index1 ] - 1. * moints[ index2 ] ) * _cvec( j ); 
         }else{

           if( pindex2 != psir ){
             if( psis != old_psis )
               read_gamma_mointb_ind_p( mointfile, psis, dble_arr2, indx_arr2, size_arr2, optimize_size_arr );
             old_psis = psis;
           }

           // now we create the first integral 
           moint1 = 0.0;
           if( pindex1 == psir ){
             INDEX_TO_BE_FOUND = index1;
             mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
             if( mo_index_in_arr >= 0 ){
               moint1 = dble_arr[ mo_index_in_arr ];
               lower_bound1 = mo_index_in_arr;
             }
           }else{
             INDEX_TO_BE_FOUND = index1;
             mo_index_in_arr = binary_search( indx_arr2, INDEX_TO_BE_FOUND, 0, size_arr2-1 );
             if( mo_index_in_arr >= 0 ){
               moint1 = dble_arr2[ mo_index_in_arr ];
               lower_bound1 = mo_index_in_arr;
             }
           }

           // now we create the second integral
           moint2 = 0.0;
           if( pindex2 == psir ){
             INDEX_TO_BE_FOUND = index2;
             mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
             if( mo_index_in_arr >= 0 ){
               moint2 = dble_arr[ mo_index_in_arr ];
               lower_bound1 = mo_index_in_arr;
             }
           }else{
             INDEX_TO_BE_FOUND = index2;
             mo_index_in_arr = binary_search( indx_arr2, INDEX_TO_BE_FOUND, 0, size_arr2-1 );
             if( mo_index_in_arr >= 0 ){
               moint2 = dble_arr2[ mo_index_in_arr ];
               lower_bound1 = mo_index_in_arr;
             }
           }


           cival += ( 2. * moint1 - 1. * moint2 ) * _cvec( j );

           old_pindex1 = pindex1;
           old_pindex2 = pindex2;
         }
      }
   sigma_vec(i) = cival;   
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
