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

int get_ijkl( int i, int j, int k, int l ){
    int pq, rs, index;
    if( i > j ) 
        pq = (i*(i+1))/2+j;
    else        
        pq = (j*(j+1))/2+i;

    if( k > l ) 
        rs = (k*(k+1))/2+l;
    else        
        rs = (l*(l+1))/2+k;

    if( pq > rs ) 
        index = (pq*(pq+1))/2+rs;
    else          
        index = (rs*(rs+1))/2+pq;

    return index;
}


void Davidson::Init(
    UnitCell& inUCell,
    SuperCell& inSCell
){
   cout << "initializing davidson..." << flush;
   full_moint_set = false;
   UCell = inUCell;
   SCell = inSCell;
   cout << "done." << endl;
}

void Davidson::set_evals_and_moints(
    VMatrixXd& eigenvals,
    vector< double >& moints_
){
    full_moint_set = true;
    evals = eigenvals;
    moints = moints_;
}

void Davidson::set_evals( VMatrixXd& eigenvals ){ evals = eigenvals; }

void Davidson::set_evals_and_moints(
    VMatrixXd& eigenvals,
    vector< complex< double > >& moints_
){
    full_moint_set = true;
    evals = eigenvals;
    kmoints = moints_;
}

int get_ijkl_kpoint_(
    int kindex,
    int orbindex,
    int max_k,
    int max_orb
){
    int index = kindex + max_k * orbindex;
    return index;
}



void Davidson::DavidsonCIS_kpoint(
    int sub_size,
    int numberOfEvals_,
    int tol_,
    bool write_eigenvalues,
    const char* filename 
){
   tol = pow( 10., 1. * tol_ );

   int psia,psib,psii,psij;
   int pq, rs, index1, index2;
   int numberOfEvals = numberOfEvals_;

   nocc = (int)UCell.nao/2;
   int nmo_ucell = UCell.nao;
   int nkpt = SCell.nkpt;
   cidim = nocc * nocc * nkpt;
   int max_k = nkpt * nkpt * nkpt;
   int max_orb = nmo_ucell * nmo_ucell * nmo_ucell * nmo_ucell;

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

   

   if( numberOfEvals >= (int)floor(cidim/2) ){
      cout << "Requested Davidson... but the number of eigenvalues needed is over half "
              "the dimension of the cis matrix... so let's just "
              "diagonalize the full matrix!" << endl;
      ofstream outstream;
      outstream.open( filename );
      for( int i = 0; i < nkpt; ++i ){
        int k0 = SCell.reduced_k[ i ]( 0 ); if( k0 > (int)floor( SCell.nkx / 2 ) ) k0 -= SCell.nkx;
        int k1 = SCell.reduced_k[ i ]( 1 ); if( k1 > (int)floor( SCell.nky / 2 ) ) k1 -= SCell.nky;
        int k2 = SCell.reduced_k[ i ]( 2 ); if( k2 > (int)floor( SCell.nkz / 2 ) ) k2 -= SCell.nkz;
        printf( "Doing k-point : %3d ... Reduced k-point (%3d,%3d,%3d) || (%3d,%3d,%3d) \n", i, 
                SCell.reduced_k[ i ]( 0 ), SCell.reduced_k[ i ]( 1 ), SCell.reduced_k[ i ]( 2 ),
                k0, k1, k2 );
        CheckerkCIS( i, write_eigenvalues, outstream );
        cout << "done." << endl;
      }
      outstream.close();
      return;
   }

   

   int cdim;
   int oldsize, currentsize, correctedsize;
   int counter = 0;
   int nmo_ucell_sq = nmo_ucell * nmo_ucell;
   int nmo_ucell_cu = nmo_ucell_sq * nmo_ucell;
   int orbindex1, kindex1, orbindex2, kindex2;
   VectorXcd   residual;   residual.resize( cidim );
   VectorXcd correction; correction.resize( cidim );
   vector< VectorXcd > HrvecsXcd;
   vector< VectorXcd > resEvecsXcd;
   VectorXcd tempEvecXcd;
   tempEvecXcd.resize( cidim );
   VectorXd  tempEvecXd;
   VectorXcd projV;  projV.resize( cidim );
   eigenpairsXcd epairs;
   HdiagXcd.resize( cidim );
   complex< double > oneXcd ( 1.0, 0.0 );
   complex< double > zeroXcd ( 0.0, 0.0 );
   MatrixXcd OldSubspaceMatr;
   MatrixXcd SubspaceMatr;
   int nconv;
   bool converged = false;
   complex< double > cival;
   int lkp;
   bool found;
   double norm;
   ofstream outfile;
   if( write_eigenvalues ){
     outfile.open( filename );
   } 

   //
   //
   // Looping over all delta k's...
   //
   //
   for( int kp = 0; kp < nkpt; ++kp ){

      // Need to clear list of 'r', 'Hr', and set number of converged eigenvalues to 0
      ckpoint = kp;
      HrvecsXcd.clear();
      resEvecsXcd.clear();
      nconv = 0;
      converged = false;

      //
      //
      // Making the diagonal matrix of H (used to make the approximate inverse of H) 
      //
      //
      //cout << "SETTING UP DIAGONAL MATRIX FOR KPOINT : '" << kp << "'..." << flush;
      printf( "-------------------------------------------\n" );
      printf( "%-50s %3d \n", "SETTING UP DIAGONAL MATRIX FOR KPOINT : ", kp );
      for( int i = 0; i < nocc * nocc; ++i ){ // Looping over the virtual and occupied space
        psii = i / nocc;
        psia = i - int( i / nocc ) * nocc + nocc;
        psij = psii;
        psib = psia;
        orbindex1 = psia + nmo_ucell * psii + nmo_ucell_sq * psij + nmo_ucell_cu * psib;
        orbindex2 = psia + nmo_ucell * psib + nmo_ucell_sq * psij + nmo_ucell_cu * psii;

        for( int ikp = 0; ikp < nkpt; ++ikp ){ // Looping over all k's
          Vector3i lkpoint = SCell.reduced_k[ ikp ] + SCell.reduced_k[ kp ];
          
          if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
          if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
          if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
          
          if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
          if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
          if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
          
          lkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );

          kindex1 = lkp + ikp * nkpt + ikp * nkpt * nkpt;
          kindex2 = lkp + lkp * nkpt + ikp * nkpt * nkpt;
          index1 = get_ijkl_kpoint_( kindex1, orbindex1, max_k, max_orb );
          index2 = get_ijkl_kpoint_( kindex2, orbindex2, max_k, max_orb );

          real( cival ) = evals.irrep( lkp )( psia, 0 ) - evals.irrep( ikp )( psii, 0 );
          imag( cival ) = 0.0;
          cival += 2. * kmoints[ index1 ] - 1. * kmoints[ index2 ]; 

          HdiagXcd( i + ikp * nocc * nocc ) = cival; 
        }  


      }

      cdim = cidim;

      //
      //
      // Let's make the random vector as a first guess ... and make it non-imaginary for everyone's sakes
      //
      //

      printf( "%-50s \n", "MAKING STARTING GUESS ... " );
      tempEvecXd = VectorXd::Random( cdim );
      //for( int i = 0; i < cdim; ++i ) tempEvecXcd( i ) = oneXcd * tempEvecXd( i );
      //tempEvecXcd = VectorXcd::Random( cdim );

      resEvecsXcd.push_back( tempEvecXcd ); 
      currentsize = 1; 
      oldsize = currentsize;
      correctedsize = currentsize;

      //
      //
      // Let's make our first vector and Hr vector if it hasn't been made already, and orthogonalizing
      // all new vectors to the space spanned by the current set of vectors.
      //
      //
      printf( "%-50s %3d \n", "ORTHOGONALIZING CURRENT SUBSPACE : SIZE = ", currentsize );
      for( int i = 0; i < currentsize; i++ ){
         if( i == 0 ){
            resEvecsXcd[ 0 ] = resEvecsXcd[ 0 ].normalized();
            HrvecsXcd.push_back( fill_hrvecXcd( resEvecsXcd[ 0 ] ) );
         }else{
            for( int j = 0; j < i; j++ ){
               resEvecsXcd[ i ] -=  ( resEvecsXcd[ j ] * (resEvecsXcd[ j ].adjoint() * resEvecsXcd[ i ] )) ;
            }
            resEvecsXcd[ i ] = resEvecsXcd[ i ].normalized();
            HrvecsXcd.push_back( fill_hrvecXcd( resEvecsXcd[ i ] ) );
         }
      }

      //
      //
      // Now we create the "little H" matrix, the hamiltonian projected onto the current basis 
      //
      //
      printf( "%-50s %3d \n", "CREATING ( b^T H b ) matrix : SIZE = ", currentsize );
      OldSubspaceMatr.resize( oldsize, oldsize );
      for( int i = 0; i < oldsize; i++ ){
         for( int j = 0; j <= i; j++ ){
            OldSubspaceMatr( i, j ) = ( resEvecsXcd[ i ].adjoint() ) * HrvecsXcd[ j ];
            OldSubspaceMatr( j, i ) = conj( OldSubspaceMatr( i, j ) );
         }
      }

      //
      //
      // Starting convergence cycle... 
      //
      //
      printf( "%-50s \n", "STARTING ITERATIONS..." );
      while( !converged ){
        if( oldsize != currentsize ){ // We again orthogonalize everything...
          printf( "%-50s %3d / %3d \n", "ORTHOGONALIZING SUBSPACE (CURRENTSIZE/OLDSIZE)", currentsize, oldsize );

          for( int i = oldsize; i < currentsize; i++ ){
             if( i == 0 ){
                resEvecsXcd[ 0 ] = resEvecsXcd[ 0 ].normalized();
                HrvecsXcd.push_back( fill_hrvecXcd( resEvecsXcd[ 0 ] ) );
             }else{
                for( int j = 0; j < i; j++ ){
                   complex< double > factor (0.0,0.0);
                   for( int k = 0; k < resEvecsXcd[ i ].size(); ++k )
                     factor += conj( resEvecsXcd[ j ]( k ) ) * resEvecsXcd[ i ]( k );
                   resEvecsXcd[ i ] -= ( (factor) * resEvecsXcd[ j ] ); 
                }
                if( resEvecsXcd[ i ].norm() < 1e-8 ){
                   cout << "CORRECTION VECTOR IS NEARLY DEGENERATE W/ CURRENT SUBSPACE (NORM = " << resEvecsXcd[ i ].norm() << ") " << endl;
                }
                while( resEvecsXcd[ i ].norm() < 1e-8 ){ // We may have numerical instabilities ... also the current vector isn't really that much 'outside' the subspace ...
                                                         // so just use a new vector!
                   cout << " o REPLACING W/ RANDOM VECTOR ... " << flush;
                   resEvecsXcd[ i ] = VectorXcd::Random( cdim );
                   resEvecsXcd[ i ] = resEvecsXcd[ i ].normalized();
                   for( int j = 0; j < i; j++ ){
                      complex< double > factor = resEvecsXcd[ j ].adjoint() * resEvecsXcd[ i ];
                      resEvecsXcd[ i ] -=  ( (factor) * resEvecsXcd[ j ] ) ;
                   }
                   cout << " NEW NORM : " << resEvecsXcd[ i ].norm() << endl;
                }
                resEvecsXcd[ i ] = resEvecsXcd[ i ].normalized();
                HrvecsXcd.push_back( fill_hrvecXcd( resEvecsXcd[ i ] ) );
             }
          }

        }

        // Checking for numerical instabilities with the Gram-Schmidt...
        double rnorm1;
        for( int i = 0; i < currentsize; ++i ){
        for( int j = 0; j < currentsize; ++j ){
          rnorm1 = ( resEvecsXcd[ i ].adjoint() * ( resEvecsXcd[ j ] ) ).norm();
          if( i == j && fabs( rnorm1 - 1.) > 1e-10 ){
            cout << "NOT NORMALIZED AFTER...( " << i << ", " << j << " ) " << endl;
            cout << rnorm1 << endl;
            exit( EXIT_FAILURE );
          }
          if( i != j && fabs( rnorm1 ) > 1e-1 ){
            cout << "NOT ORTHOGONAL AFTER...( " << i << ", " << j << " ) " << endl;
            cout << rnorm1 << endl;
            exit( EXIT_FAILURE );
          }
        }
        }

        cout << "CREATING LITTLE H [size= " << currentsize << "]... " << flush;
        SubspaceMatr.resize( currentsize, currentsize );
        SubspaceMatr.block( 0, 0, oldsize, oldsize ) = OldSubspaceMatr; // Reuse parts of old subspace matrix
        for( int i = oldsize; i < currentsize; i++ ){
           for( int j = 0; j <= i; j++ ){
              SubspaceMatr( i, j ) = ( resEvecsXcd[ i ].adjoint() ) * HrvecsXcd[ j ];
              SubspaceMatr( j, i ) = conj( SubspaceMatr(i,j) );
           }
        }
        OldSubspaceMatr.resize( currentsize, currentsize );
        OldSubspaceMatr = SubspaceMatr;
        oldsize = currentsize;
        cout << "done." << endl;

        //
        //
        // Let's find the eigenvalues and vectors for this subspace.
        //
        //
        cout << "FINDING EIGENVALUES IN SUBSPACE ... " << flush;
        SelfAdjointEigenSolver< MatrixXcd > eigensolver( SubspaceMatr );
        epairs.evecs = eigensolver.eigenvectors();
        epairs.evals = eigensolver.eigenvalues();
        cout << "done." << endl;

        /*
        cout << "CHECKING IF SOLUTION EXISTS..." << endl;
        bool solution_exists = ( epairs.evecs * epairs.evals.asDiagonal() * epairs.evecs.adjoint() ).isApprox( SubspaceMatr, 1e-10 );
        cout << "THE 'SOLUTION' : " << endl;
        cout << epairs.evecs * epairs.evals.asDiagonal() * epairs.evecs.adjoint() << endl;
        cout << "SUBSPACE MATR : " << endl;
        cout << SubspaceMatr << endl;
        if( !solution_exists ){
          cout << "BAD SOLUTION ..." << endl;
          exit( EXIT_FAILURE );
        }
        */

        
        //
        // 
        // Now calculating the residual vector, finding the eigenfunction in the original basis
        // and calculating |A*v - e*v|.  Assuming we have already converged the lowest (m-1) eigenvalues and eigenvectors, 
        // we now pick out the m-th eigenvalue and eigenvector 
        //
        //

        cout << "CALCULATING RESIDUAL ... " << flush;
        residual *= 0.0;
        for( int i = 0; i < cdim; i++ ){
           for( int j = 0; j < currentsize; j++ ){
              residual( i ) += HrvecsXcd[ j ]( i ) * epairs.evecs( j, nconv ); // H * v
              residual( i ) -= epairs.evals( nconv ) * resEvecsXcd[ j ]( i ) * epairs.evecs( j, nconv ); // - e * r
           }
        }
        cout << "res = " << residual.norm() << endl;

        if( ( residual.norm() ) > tol){ // Checking to see whether the residual is small enough
           if( currentsize == (sub_size-1) ){ // Checks if we can increase the size of our subspace
              printf("Could not find lowest [%d] eigenvalues with MaxSubspaceSize[%d]\n",numberOfEvals,sub_size);
              converged = 1;
           }
           else if( currentsize == cidim ){ // we reached the maximum subspace... so we definitely converged!
              cout << "Residual not reached... but reached size of CI matrix, all eigenvalues are converged!..." << endl;
              cout << "lowest eigenvalues (H, eV):" << endl;
              for ( int i = 0; i < numberOfEvals; i++ ){
                 printf("%8.10e, %8.10e \n",  epairs.evals( i ), ( epairs.evals( i ) * 27.211 ) );
                 if( write_eigenvalues ){ 
                   outfile << fixed << setprecision( 16 ) << epairs.evals( i ) << endl;
                 }
              }
              cout << "converged with dimension : " << currentsize << endl;
              fflush(stdout); 
              converged = true;
           }
           else{ // We compute the correction vector and add it to our subspace 
              if( residual.norm() < 1e-6 ){
                for( int i = 0; i < cdim; i++ ){
                   correction( i ) = ( - 1. / ( real( HdiagXcd( i ) ) - epairs.evals( nconv ) ) ) * residual( i );
                }
              }else{
                double lowest_eigen = 1e-6;
                for( int i = nconv; i < currentsize; ++i ){
                  if( epairs.evals( i ) > lowest_eigen ){
                    lowest_eigen = epairs.evals( i );
                    break;
                  } 
                }
                for( int i = 0; i < cdim; i++ ){
                   correction( i ) = ( - 1. / ( real( HdiagXcd( i ) ) - lowest_eigen ) ) * residual( i );
                }
              }
              resEvecsXcd.push_back( correction );
              correctedsize++;
              printf("Current eval [%d] : %20.16e \n", nconv, epairs.evals( nconv ));
              fflush(stdout);
           }
        }
       
        //If converged, see whether to go on to a new eigenvalue if it was specified
        if( (residual.norm() ) < tol){
           if( nconv == ( numberOfEvals - 1 ) ){
              cout << "lowest eigenvalues (H, eV):" << endl;
              for ( int i = 0; i < numberOfEvals; i++ ){
                 printf("%20.16f, %20.16f \n",  epairs.evals( i ), ( epairs.evals( i ) * 27.211 ) );
                 if( write_eigenvalues ){
                   outfile << fixed << setprecision( 16 ) << epairs.evals( i ) << endl; 
                 }
              }
              cout << "converged with dimension : " << currentsize << endl;
              fflush(stdout);
              converged = true;
           }
           else if( currentsize == cidim ){ // we reached the maximum subspace... so we definitely converged!
              cout << "Residual reached.. and reached size of CI matrix, all eigenvalues are converged!..." << endl;
              cout << "lowest eigenvalues (H, eV):" << endl;
              for ( int i = 0; i < numberOfEvals; i++ ){
                 printf("%20.16f, %20.16f \n",  epairs.evals( i ), ( epairs.evals( i ) * 27.211 ) );
                 if( write_eigenvalues ){ 
                   outfile << fixed << setprecision( 16 ) << epairs.evals( i ) << endl;
                 }
              }
              cout << "converged with dimension : " << currentsize << endl;
              fflush(stdout); 
              converged = true;
           }
           else{
              printf("RESIDUAL < TOL : CURRENT EVAL [%d] CONVERGED => %14.6e \n", nconv, epairs.evals( nconv ) );
              nconv++;
              // 
              // 
              // If our current subspace of size N finds the lowest N eigenvalues up to desired accuracy,
              // need to increase the size of our space done so by adding in a random vector
              // 
              // 
              if( nconv >= currentsize ){
                tempEvecXd = VectorXd::Random( cdim );
                for( int i = 0; i < cdim; ++i ) tempEvecXcd( i ) = oneXcd * tempEvecXd( i );
                resEvecsXcd.push_back( tempEvecXcd );
                correctedsize++;
              }
           }
        }
       
        currentsize = correctedsize;   
      } // end convergence loop

   } // end loop over delta k's
   outfile.close();

}



VectorXcd Davidson::fill_hrvecXcd( VectorXcd& cvec ){
   VectorXcd sigma_vec;
   sigma_vec.resize( cidim );

   int psii, psij, psia, psib;
   int orbindex1, orbindex2;
   int kindex1, kindex2, index1, index2;
   int nmo_ucell = UCell.nao;
   int nmo_ucell_sq = nmo_ucell * nmo_ucell;
   int nmo_ucell_cu = nmo_ucell * nmo_ucell_sq;
   int nocc = (int)(UCell.nao/2);
   int nkpt = SCell.nkpt;
   int max_k = nkpt * nkpt * nkpt;
   int max_orb = nmo_ucell_cu * nmo_ucell; 
   int noccsq = nocc * nocc;
   int sindex1, sindex2;
   complex< double > cival;
   complex< double > oneXcd ( 1.0, 0.0 );
   int lkp, mkp;
   bool found;

   for( int i = 0; i < noccsq; ++i ){ // Looping over the virtual and occupied space for left side
     psii = i / nocc;
     psia = i - int( i / nocc ) * nocc + nocc;
     for( int ikp = 0; ikp < nkpt; ++ikp ){ // Looping over all k1's 

       sindex1 = i + ikp * noccsq;

       Vector3i lkpoint = SCell.reduced_k[ ikp ] + SCell.reduced_k[ ckpoint ];

       if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
       if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
       if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
  
       if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
       if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
       if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
  
       lkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );

       cival = oneXcd * ( evals.irrep( lkp )( psia, 0 ) - evals.irrep( ikp )( psii, 0 ) ) * cvec( sindex1 );
//       cout << "MAKING SIGMA I,A (kindex = " << ikp << " ): " << psii << ", " << psia << endl;
       for( int j = 0; j < noccsq; ++j ){ // Looping over the virtual and occupied space for right side
         psij = j / nocc;
         psib = j - int( j / nocc ) * nocc + nocc;
         orbindex1 = psia + nmo_ucell * psii + nmo_ucell_sq * psij + nmo_ucell_cu * psib;
         orbindex2 = psia + nmo_ucell * psib + nmo_ucell_sq * psij + nmo_ucell_cu * psii;
      
         for( int jkp = 0; jkp < nkpt; ++jkp ){ // Looping over all k2's 
           lkpoint = SCell.reduced_k[ jkp ] + SCell.reduced_k[ ckpoint ];
          
           if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
           if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
           if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
          
           if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
           if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
           if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
          
           mkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );
      
           kindex1 = lkp + ikp * nkpt + jkp * nkpt * nkpt;
           kindex2 = lkp + mkp * nkpt + jkp * nkpt * nkpt;
           index1 = get_ijkl_kpoint_( kindex1, orbindex1, max_k, max_orb );
           index2 = get_ijkl_kpoint_( kindex2, orbindex2, max_k, max_orb );
      
           sindex2 = j + jkp * noccsq;
      
//           cout << "ADDING A,I,J,B k1,k2 : " << psia << ", " << psii << ", " << psij << ", " << psib << ", " << ikp << ", " << jkp << endl;
           cival += ( 2. * kmoints[ index1 ] - 1. * kmoints[ index2 ] ) * cvec( sindex2 ); 
         }
       }
       sigma_vec( sindex1 ) = cival;
        
     }  
   }

   return sigma_vec;
}



void Davidson::CheckerkCIS( 
   int inkpoint,
   bool write_eigenvalues,
   std::ofstream& outstream
 ){
   int psii, psij, psia, psib;
   int orbindex1, orbindex2;
   int kindex1, kindex2, index1, index2;
   int nmo_ucell = UCell.nao;
   int nmo_ucell_sq = nmo_ucell * nmo_ucell;
   int nmo_ucell_cu = nmo_ucell * nmo_ucell_sq;
   int nocc = (int)(UCell.nao/2);
   int nkpt = SCell.nkpt;
   int max_k = nkpt * nkpt * nkpt;
   int max_orb = nmo_ucell_cu * nmo_ucell; 
   int noccsq = nocc * nocc;
   int sindex1, sindex2;
   complex< double > cival;
   complex< double > oneXcd ( 1.0, 0.0 );
   complex< double > zeroXcd ( 0.0, 0.0 );
   int lkp, mkp;
   bool found;
   MatrixXcd cimatr;
   cimatr.resize( noccsq * nkpt, noccsq * nkpt );
   ckpoint = inkpoint;
   for( int j = 0; j < cimatr.cols(); ++j )
   for( int i = 0; i < cimatr.rows(); ++i )
     cimatr( i, j ) = zeroXcd;

   for( int i = 0; i < noccsq; ++i ){ // Looping over the virtual and occupied space for left side
     psii = i / nocc;
     psia = i - int( i / nocc ) * nocc + nocc;
     for( int ikp = 0; ikp < nkpt; ++ikp ){ // Looping over all k1's 

       sindex1 = i + ikp * noccsq;

       Vector3i lkpoint = SCell.reduced_k[ ikp ] + SCell.reduced_k[ ckpoint ];

       if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
       if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
       if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
  
       if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
       if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
       if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
  
       lkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );

       cimatr( sindex1, sindex1 ) = oneXcd * ( evals.irrep( lkp )( psia, 0 ) - evals.irrep( ikp )( psii, 0 ) );
//       cout << "MAKING SIGMA I,A (kindex = " << ikp << " ): " << psii << ", " << psia << endl;
       for( int j = 0; j < noccsq; ++j ){ // Looping over the virtual and occupied space for right side
         psij = j / nocc;
         psib = j - int( j / nocc ) * nocc + nocc;
         orbindex1 = psia + nmo_ucell * psii + nmo_ucell_sq * psij + nmo_ucell_cu * psib;
         orbindex2 = psia + nmo_ucell * psib + nmo_ucell_sq * psij + nmo_ucell_cu * psii;
      
         for( int jkp = 0; jkp < nkpt; ++jkp ){ // Looping over all k2's 
           lkpoint = SCell.reduced_k[ jkp ] + SCell.reduced_k[ ckpoint ];
          
           if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
           if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
           if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
          
           if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
           if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
           if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
          
           mkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );
      
           kindex1 = lkp + ikp * nkpt + jkp * nkpt * nkpt;
           kindex2 = lkp + mkp * nkpt + jkp * nkpt * nkpt;
           index1 = get_ijkl_kpoint_( kindex1, orbindex1, max_k, max_orb );
           index2 = get_ijkl_kpoint_( kindex2, orbindex2, max_k, max_orb );
      
           sindex2 = j + jkp * noccsq;
      
           cimatr( sindex1, sindex2 ) += ( 2. * kmoints[ index1 ] - 1. * kmoints[ index2 ] );
         }
       }
     }  
   }

   SelfAdjointEigenSolver< MatrixXcd > solver( cimatr );
   cout << "EIGENVALUES : " << endl;
   for( int i = 0; i < cidim; ++i ) printf( "%20.16f \n", solver.eigenvalues()( i ) );

   if( write_eigenvalues ){
     for( int i = 0; i < noccsq * nkpt; ++i )
       outstream << fixed << setprecision( 16 ) << solver.eigenvalues()( i ) << endl;
   }

/*
   cout << "EIGENVECTOR 0 : " << endl;
   cout << solver.eigenvectors().col( 0 ) << endl;
   if( cimatr.rows() > 1 ){
     cout << "EIGENVECTOR 1 : " << endl;
     cout << solver.eigenvectors().col( 1 ) << endl;
   }
*/
}


