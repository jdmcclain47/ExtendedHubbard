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
#include "mpi.h"
#include "binary_search.h"
#include "moint_read_file.h"
#include "mointegrals.h"
#include "common.h"
#include "sorter.h"
#include "cellinfo.h"

using namespace std;
using namespace Eigen;

static const int WORK_TAG = 1;
static const int DIE_TAG = 2;
static bool is_initialized = false;
static Eigen::MatrixXd corr1matr;
static Eigen::MatrixXd corr2matr;
static Eigen::MatrixXd evecsXd;
static std::vector< double > evals;

static size_t get_ijkl( int i, int j, int k, int l ){
    size_t pq, rs, index;
    pq = ( i > j ) ? (i*(i+1))/2+j : (j*(j+1))/2+i;
    rs = ( k > l ) ? (k*(k+1))/2+l : (l*(l+1))/2+k;
    index = ( pq > rs ) ? (pq*(pq+1))/2+rs : (rs*(rs+1))/2+pq; 
    return index;
}

static void make_corr_matr( 
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints
){
    int nmo_ucell = UCell.nao;
    int ntrans = SCell.total_number_of_cells;
    int index;
    bool found;

    for( int iao = 0; iao < nmo_ucell; ++iao ){
    for( int jao = 0; jao < nmo_ucell; ++jao ){
      for( int itrans = 0; itrans < ntrans; ++itrans ){
      for( int jtrans = 0; jtrans < ntrans; ++jtrans ){
        Vector3i tvec = SCell.reduced_t[ jtrans ] - SCell.reduced_t[ itrans ];
        found = false;
        index = SCell.get_supercell_index( tvec( 0 ), tvec( 1 ), tvec( 2 ), &found );

        corr1matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) = aoints.getCorr1Int( iao, jao, index );
        corr2matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) = aoints.getCorr2Int( iao, jao, index );
        if( iao == jao && itrans == jtrans ){ // we are onsite
          corr1matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getHubbardU( SCell.ElementStr[ iao ] );
          corr2matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getHubbardU( SCell.ElementStr[ iao ] );
        }
      }
      }
    }
    }
}


  
double get_moint( 
    vector< double >& dble_arr,
    vector< size_t >& indx_arr,
    size_t size_arr,
    int p, int q, int r, int s
){
    double outval = 0.0;
    size_t INDEX_TO_BE_FOUND = get_ijkl( p, q, r, s );
    int mo_index_in_arr = binary_search( indx_arr, INDEX_TO_BE_FOUND, 0, size_arr-1 );
    outval = ( mo_index_in_arr >= 0 ) ? dble_arr[ mo_index_in_arr ] : 0.0;
    return outval;
}


void create_hrvec(
    VectorXd& hrvec,
    VectorXd& rvec,
    const int nocc,
    const int norb
){
  /* Setting everything up for MOINTEGRALS */
  vector< double > dble_arr1, dble_arr2;
  vector< size_t > indx_arr1, indx_arr2;
  size_t nmax = 1000000;
  size_t arr_size1, arr_size2;
  size_t opt_size_arr = 0;
  bool optimize_size_arr = true;

  indx_arr1.resize( nmax );
  dble_arr1.resize( nmax );
  indx_arr2.resize( nmax );
  dble_arr2.resize( nmax );
  /* DONE setting everything up for MOINTEGRALS */

  MPI_Status status;
  vector< double > t_arr;
  int psia, psib, psii, psij, ia_index, jb_index;
  int taskid, numtasks, nwork, currentwork, which_work;
  double moint1, moint2;
  bool done = false;
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
  nwork = norb - nocc;
  //
  //
  // SETTING UP MASTER
  //
  //
  if( taskid == 0 ){
    //
    //
    // sending out initial workloads to the slaves
    //
    //
    currentwork = 0;
    for( int rank = 1; rank < numtasks; ++rank ){
      //cout << "SENDING OUT SOME WORK TO RANK " << rank << " :) !!! work = " << currentwork << endl;
      MPI_Send( &currentwork, 1, MPI_INT, rank, WORK_TAG, MPI_COMM_WORLD );
      currentwork++;
    }
    //
    //
    // receive results and send out more work if needed
    //
    //
    while( currentwork < nwork ){
      //cout << "TASK 0 : RECEIVING DATA...." << endl;
      MPI_Recv( &which_work, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      MPI_Recv( &hrvec( which_work * nocc ), nocc, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      //cout << "   CAME FROM " << status.MPI_SOURCE << endl;
      //cout << "   SENDING OUT WORK =  " << currentwork << " TO TASK : " << status.MPI_SOURCE << endl;
      MPI_Send( &currentwork, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD ); 
      currentwork++;
    }
    //
    //
    // receive any outstanding results
    //
    //
    for( int rank = 1; rank < numtasks; ++rank ){
      MPI_Recv( &which_work, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      //cout << "TASK 0 RECIEVED WORK = " << which_work << endl;
      if( which_work < nwork )
        MPI_Recv( &hrvec( which_work * nocc ), nocc, MPI_DOUBLE, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      MPI_Send( 0, 0, MPI_INT, rank, DIE_TAG, MPI_COMM_WORLD );
    }
    //cout << "HRVEC " << endl; cout << hrvec << endl;
  }
  //
  //
  // SETTING UP SLAVES 
  //
  //
  else{
    t_arr.resize( nocc );
    while( !done ){
      //cout << "TASK " << taskid << " RECEIVING SOME WORK M8!! " << flush; 
      MPI_Recv( &which_work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
      //cout << "WORK = " << which_work << endl; 
      if( status.MPI_TAG == DIE_TAG /*|| which_work == nwork*/ ){
        //cout << "TASK " << taskid << " JK M8, DYING....!! " << endl; 
        break;
      }
      //
      //
      // creating the sigma_ia = eval(a)-eval(i)+sum_{jb}c_bj*[2*(ai|jb)-(ab|ji)]
      //
      //
      psia = which_work + nocc;
      gamma_moint_to_vector( dble_arr1, indx_arr1, corr1matr, evecsXd, norb,
         arr_size1, psia, (psia+1), 0, nocc, 0, nocc, nocc, norb, -10 ); 
      //gamma_moint_to_vector( dble_arr1, indx_arr1, corr1matr, evecsXd, norb,
      //   arr_size1, nocc, norb, 0, nocc, 0, nocc, nocc, norb, -10 ); 
      sorter( dble_arr1, indx_arr1, arr_size1 );
      gamma_moint_to_vector( dble_arr2, indx_arr2, corr2matr, evecsXd, norb,
         arr_size2, psia, (psia+1), nocc, norb, 0, nocc, 0, nocc, -10 ); 
      //gamma_moint_to_vector( dble_arr2, indx_arr2, corr2matr, evecsXd, norb,
      //   arr_size2, nocc, norb, nocc, norb, 0, nocc, 0, nocc, -10 ); 
      sorter( dble_arr2, indx_arr2, arr_size2 );
      for( psii = 0; psii < nocc; ++psii ){
      ia_index = ( psia - nocc ) * nocc + psii;
      t_arr[ psii ] = ( evals[ psia ] - evals[ psii ] ) * rvec( ia_index );
        //cout << "FOR PSIA, PSII : " << psia << ", " << psii << endl;
      for( psij = 0; psij < nocc; ++psij ){
      for( psib = nocc; psib < norb; ++psib ){
        jb_index = ( psib - nocc ) * nocc + psij;
        moint1 = get_moint( dble_arr1, indx_arr1, arr_size1, psia, psii, psij, psib );
        moint2 = get_moint( dble_arr2, indx_arr2, arr_size2, psia, psib, psij, psii );
        t_arr[ psii ] += ( 2. * moint1 - moint2 ) * rvec( jb_index ); 
        //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psii, psij, psib, moint1 );
        //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psib, psij, psii, moint2 );
      }
      }
      }
      //cout << "TASK " << taskid << " DONE WITH WORK, SENDING BACK TO 0!! " << endl; 
      MPI_Send( &which_work, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ); 
      MPI_Send( &t_arr[ 0 ], nocc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD ); 
    } 
  }

}


void create_hdiag( 
    VectorXd& Hdiag,
    const int nocc,
    const int norb
){
  /* Setting everything up for MOINTEGRALS */
  vector< double > dble_arr1, dble_arr2;
  vector< size_t > indx_arr1, indx_arr2;
  size_t nmax = 1000000;
  size_t arr_size1, arr_size2;
  size_t opt_size_arr = 0;
  bool optimize_size_arr = true;

  indx_arr1.resize( nmax );
  dble_arr1.resize( nmax );
  indx_arr2.resize( nmax );
  dble_arr2.resize( nmax );
  /* DONE setting everything up for MOINTEGRALS */

  MPI_Status status;
  vector< double > t_arr;
  int psia, psib, psii, psij;
  int taskid, numtasks, nwork, currentwork, which_work;
  double moint1, moint2;
  bool done = false;
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
  nwork = norb - nocc;
  //
  //
  // SETTING UP MASTER
  //
  //
  if( taskid == 0 ){
    //
    //
    // sending out initial workloads to the slaves
    //
    //
    currentwork = 0;
    for( int rank = 1; rank < numtasks; ++rank ){
      //cout << "SENDING OUT SOME WORK TO RANK " << rank << " :) !!! work = " << currentwork << endl;
      MPI_Send( &currentwork, 1, MPI_INT, rank, WORK_TAG, MPI_COMM_WORLD );
      currentwork++;
    }
    //
    //
    // receive results and send out more work if needed
    //
    //
    while( currentwork < nwork ){
      //cout << "TASK 0 : RECEIVING DATA...." << endl;
      MPI_Recv( &which_work, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      MPI_Recv( &Hdiag( which_work * nocc ), nocc, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      //cout << "   CAME FROM " << status.MPI_SOURCE << endl;
      //cout << "   SENDING OUT WORK =  " << currentwork << " TO TASK : " << status.MPI_SOURCE << endl;
      MPI_Send( &currentwork, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD ); 
      currentwork++;
    }
    //
    //
    // receive any outstanding results
    //
    //
    for( int rank = 1; rank < numtasks; ++rank ){
      MPI_Recv( &which_work, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      //cout << "TASK 0 RECIEVED WORK = " << which_work << endl;
      if( which_work < nwork )
        MPI_Recv( &Hdiag( which_work * nocc ), nocc, MPI_DOUBLE, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      MPI_Send( 0, 0, MPI_INT, rank, DIE_TAG, MPI_COMM_WORLD );
    }
  }
  //
  //
  // SETTING UP SLAVES 
  //
  //
  else{
    t_arr.resize( nocc );
    while( !done ){
      //cout << "TASK " << taskid << " RECEIVING SOME WORK M8!! " << flush; 
      MPI_Recv( &which_work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
      //cout << "WORK = " << which_work << endl; 
      if( status.MPI_TAG == DIE_TAG /*|| which_work == nwork*/ ){
        //cout << "TASK " << taskid << " JK M8, DYING....!! " << endl; 
        break;
      }
      //
      //
      // creating the actual diagonal part
      //
      //
      psia = which_work + nocc;
      gamma_moint_to_vector( dble_arr1, indx_arr1, corr1matr, evecsXd, norb,
         arr_size1, psia, (psia+1), 0, nocc, psia, (psia+1), 0, nocc, -10 ); 
      sorter( dble_arr1, indx_arr1, arr_size1 );
      gamma_moint_to_vector( dble_arr2, indx_arr2, corr2matr, evecsXd, norb,
         arr_size2, psia, (psia+1), psia, (psia+1), 0, nocc, 0, nocc, -10 ); 
      sorter( dble_arr2, indx_arr2, arr_size2 );
      for( psii = 0; psii < nocc; ++psii ){
        moint1 = get_moint( dble_arr1, indx_arr1, arr_size1, psia, psii, psii, psia );
        moint2 = get_moint( dble_arr2, indx_arr2, arr_size2, psia, psia, psii, psii );
        t_arr[ psii ] = evals[ psia ] - evals[ psii ];
        t_arr[ psii ] += 2. * moint1 - moint2; 
        t_arr[ psii ] = 0.0;
        //cout << "FOR PSIA, PSII : " << psia << ", " << psii << endl;
        //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psii, psii, psia, moint1 );
        //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psia, psii, psii, moint2 );
      }
      //cout << "TASK " << taskid << " DONE WITH WORK, SENDING BACK TO 0!! " << endl; 
      MPI_Send( &which_work, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ); 
      MPI_Send( &t_arr[ 0 ], nocc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD ); 
    } 
  }
        
}

void orth_vecs(
    vector< VectorXd >& vecs,
    int old_size,
    int new_size
){
    for( int i = old_size; i < new_size; ++i ){
      for( int j = 0; j < i; j++ ){
         vecs[i] -=  (vecs[j] * vecs[i].transpose() * vecs[j] ) ;
      }
      vecs[i] = vecs[i].normalized();
    }
}

void set_up_globals( 
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints, 
    int nocc,
    int norb,
    VMatrixXd& inevecs,
    VMatrixXd& inevals
){
    int numtasks, taskid;
    MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
    MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
    if( !is_initialized ){
      corr1matr.resize( norb, norb );
      corr2matr.resize( norb, norb );
      evecsXd.resize( norb, norb );
      evals.resize( norb );
      if( taskid == 0 ){
        make_corr_matr( UCell, SCell, aoints );
        evecsXd = inevecs.irrep( 0 );
        for( int i = 0; i < norb; ++i ) evals[ i ] = inevals.irrep( 0 )( i, 0 );
      }
      MPI_Bcast( &corr1matr(0,0), norb*norb, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( &corr2matr(0,0), norb*norb, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( &evecsXd(0,0), norb*norb, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( &evals[0], norb, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }
}

/* Can make the following a little bit better for when finding multiple eigenvalues....
 * After converging to an eigenvalue, add the "true" rvec to your list of rvecs,
 * then start over starting with 'nconverged' eigenvectors/eigenvalues
 */

void Davidson_Gamma_MPI(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXd& inevecs,
    VMatrixXd& inevals,
    const int numberOfEvals,
    const int itol
){

  /* Setting everything needed by all procs */
  int numtasks, taskid;
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  assert( numtasks > 1 && "DAVIDSON_GAMMA_MPI REQUIRES MORE THAN ONE PROC!" );
  int nocc = (int)(SCell.nao/2);
  int cidim = nocc * nocc;
  int norb = SCell.nao;
  MPI_Bcast( &nocc, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &norb, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &cidim, 1, MPI_INT, 0, MPI_COMM_WORLD );
  set_up_globals( UCell, SCell, aoints, nocc, norb, inevecs, inevals );

  /* Setting everything up for MOINTEGRALS */
  vector< double > dble_arr1, dble_arr2;
  vector< size_t > indx_arr1, indx_arr2;
  size_t INDEX_TO_BE_FOUND, lower_bound1, lower_bound2;
  int mo_index_in_arr;
  size_t nmax = 1000000;
  size_t opt_size_arr = 0;
  bool optimize_size_arr = true;

  indx_arr1.resize( nmax );
  dble_arr1.resize( nmax );
  indx_arr2.resize( nmax );
  dble_arr2.resize( nmax );

  /* Setting everything up for this function */
  int psia,psib,psii,psij;
  double moint1, moint2, cival;
  int pq, rs, index1, index2;
  int neval;
  double dtol = pow( 10., -1. * itol );
  clock_t t;
  int index_ci;
  VectorXd Hdiag;

  neval = numberOfEvals;
  if( numberOfEvals > cidim ){
     cout << "***********************************************************************" << endl;
     cout << "Trying to find more eigenvalues than what's possible in this subspace.." << endl;
     cout << "  Will instead find max number of eigenvalues.                         " << endl;
     cout << "***********************************************************************" << endl;
     neval = cidim;
  }
   
   // Making the diagonal Matrix
   MPI_Barrier( MPI_COMM_WORLD );
   std::cout << "CREATING DIAGONAL MATRIX..." << std::endl;
   Hdiag.resize(cidim);
   create_hdiag( Hdiag, nocc, norb );
   std::cout << "DONE CREATING DIAGONAL MATRIX..." << std::endl;

   int old_size, new_size, isconverged, error, nconverged;
   int step, max_step, n_added_states;
   double res;
   vector< VectorXd > rvec, Hrvec;
   vector< double > converged_evals;
   VectorXd tempvec, resvec, corrvec;
   MatrixXd old_sub_matr;
   MatrixXd new_sub_matr;
   eigenpairs epairs;
   tempvec.resize( cidim );
   resvec.resize( cidim );
   corrvec.resize( cidim );
   old_size = 0;
   new_size = 1; 

   max_step = 400;
   isconverged = 0;
   error = 0;
   step = 0;
   nconverged = 0;
   while( (!isconverged) && (!error) && step < max_step ){
     for( int i = old_size; i < new_size; ++i ){
       if( step == 0 )
         tempvec = VectorXd::Random( cidim );
       else
         tempvec = corrvec;
       tempvec = tempvec.normalized();
       rvec.push_back( tempvec );
     }
     orth_vecs( rvec, old_size, new_size );
     for( int i = old_size; i < new_size; ++i ){
       if( step == 0 )
         MPI_Bcast( &rvec[ i ]( 0 ), cidim, MPI_DOUBLE, 0, MPI_COMM_WORLD );
       Hrvec.push_back( rvec[ i ] );
       create_hrvec( Hrvec[ i ], rvec[ i ], nocc, norb );
       //if( taskid == 0 ){
       //  cout << "RVEC" << endl;
       //  cout << rvec[ i ].transpose() << endl;
       //  cout << "HRVEC" << endl;
       //  cout << Hrvec[ i ].transpose() << endl;
       //}
     }
     if( taskid == 0 ){
       new_sub_matr.resize( new_size, new_size );
       if( step > 0 )
         new_sub_matr.block( 0, 0, old_size, old_size ) = old_sub_matr;
       for( int j = 0; j < new_size; ++j ){
         for( int i = 0; i < j; ++i ){
           new_sub_matr( i, j ) = ( rvec[ i ].transpose() ) * Hrvec[ j ];
           new_sub_matr( j, i ) = new_sub_matr( i, j ); 
         }
         new_sub_matr( j, j ) = ( rvec[ j ].transpose() ) * Hrvec[ j ];
       }
       old_sub_matr.resize( new_size, new_size );
       old_sub_matr = new_sub_matr;

       SelfAdjointEigenSolver<MatrixXd> eigensolver( new_sub_matr );
       epairs.evals = eigensolver.eigenvalues();
       epairs.evecs = eigensolver.eigenvectors();

       for( int i = 0; i < cidim; i++){
       resvec( i ) = 0.0;
       for( int j = 0; j < new_size; j++){
         resvec(i) += Hrvec[ j ]( i ) * epairs.evecs( j, min( (new_size - 1), nconverged) );
         resvec(i) -= epairs.evals( min( (new_size - 1), nconverged) ) * rvec[j](i) * epairs.evecs( j, min( (new_size - 1), nconverged) );
       }
       }
     }

     old_size = new_size;
     new_size = new_size+1;

     if( taskid == 0 ){
       res = resvec.blueNorm();
       if( step == 0 )
         printf( "%-10s %-10s %20s \t %20s \n", "STEP", "ROOT", "EIGENVALUE", "RES" );
       printf( "  %-10d %-10d %20.14e \t %20.14e \n", (step+1), (nconverged+1), epairs.evals( nconverged ), res );
       if( res > dtol ){
         for( int i = 0; i < cidim; ++i ){
           corrvec( i ) = -1.0/( Hdiag( i ) - epairs.evals( nconverged ) ) * resvec( i );
         }  
         MPI_Bcast( &corrvec( 0 ), cidim, MPI_DOUBLE, 0, MPI_COMM_WORLD );
       }
       else{
         converged_evals.push_back( epairs.evals( nconverged ) );
         corrvec = VectorXd::Random( cidim );
         MPI_Bcast( &corrvec( 0 ), cidim, MPI_DOUBLE, 0, MPI_COMM_WORLD );
         nconverged++;
       }
       if( nconverged == cidim ) isconverged = 1;
       if( nconverged == neval ) isconverged = 1;
       if( new_size == (cidim+1) ) isconverged = 1;
     }

     if( taskid != 0 )
       MPI_Bcast( &corrvec( 0 ), cidim, MPI_DOUBLE, 0, MPI_COMM_WORLD );
     MPI_Bcast( &isconverged, 1, MPI_INT, 0, MPI_COMM_WORLD );

     step++;
   }
   if( taskid == 0 ){
     printf( "%-20s %10s %20s \n", "===============", "CIS LOWEST ROOTS", "===============" );
     printf( "  %-10s %20s \n", "ROOT", "EIGENVALUE" );
     for( int i = 0; i < nconverged; ++i ) 
       printf( "  %-10d %20.14e \n", (i+1), converged_evals[ i ] );
     printf( "%-20s %10s %20s \n", "===============", "================", "===============" );
   }
   
}


void MP2_Gamma_MPI(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXd& inevecs,
    VMatrixXd& inevals,
    const int itol
){

  /* Setting everything needed by all procs */
  MPI_Status status;
  int numtasks, taskid;
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  assert( numtasks > 1 && "MP2_GAMMA_MPI REQUIRES MORE THAN ONE PROC!" );
  int nocc = (int)(SCell.nao/2);
  int norb = SCell.nao;
  MPI_Bcast( &nocc, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &norb, 1, MPI_INT, 0, MPI_COMM_WORLD );
  set_up_globals( UCell, SCell, aoints, nocc, norb, inevecs, inevals );

  /* Setting everything up for MOINTEGRALS */
  vector< double > dble_arr1, dble_arr2;
  vector< size_t > indx_arr1, indx_arr2;
  size_t INDEX_TO_BE_FOUND, lower_bound1, lower_bound2, arr_size1, arr_size2;
  int mo_index_in_arr;
  size_t nmax = 1000000;
  size_t opt_size_arr = 0;
  bool optimize_size_arr = true;

  indx_arr1.resize( nmax );
  dble_arr1.resize( nmax );
  indx_arr2.resize( nmax );
  dble_arr2.resize( nmax );

  /* Setting everything up for this function */
  bool done = false;
  int psia,psib,psii,psij;
  int currentwork, which_work;
  double moint1, moint2, denom;
  double temp_sum = 0.0;
  double total_sum = 0.0;
  int pq, rs, index1, index2;
  int nwork = norb - nocc;
  //
  //
  // SETTING UP MASTER
  //
  //
  if( taskid == 0 ){
    //
    //
    // sending out initial workloads to the slaves
    //
    //
    currentwork = 0;
    for( int rank = 1; rank < numtasks; ++rank ){
      //cout << "SENDING OUT SOME WORK TO RANK " << rank << " :) !!! work = " << currentwork << endl;
      MPI_Send( &currentwork, 1, MPI_INT, rank, WORK_TAG, MPI_COMM_WORLD );
      currentwork++;
    }
    //
    //
    // receive results and send out more work if needed
    //
    //
    while( currentwork < nwork ){
      //cout << "TASK 0 : RECEIVING DATA...." << endl;
      MPI_Recv( &which_work, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      MPI_Recv( &temp_sum, 1, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      total_sum += temp_sum;
      //cout << "   CAME FROM " << status.MPI_SOURCE << endl;
      //cout << "   SENDING OUT WORK =  " << currentwork << " TO TASK : " << status.MPI_SOURCE << endl;
      MPI_Send( &currentwork, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD ); 
      currentwork++;
    }
    //
    //
    // receive any outstanding results
    //
    //
    for( int rank = 1; rank < numtasks; ++rank ){
      MPI_Recv( &which_work, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
      //cout << "TASK 0 RECIEVED WORK = " << which_work << endl;
      if( which_work < nwork ){
        MPI_Recv( &temp_sum, 1, MPI_DOUBLE, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
        total_sum += temp_sum;
      }
      MPI_Send( 0, 0, MPI_INT, rank, DIE_TAG, MPI_COMM_WORLD );
    }
    printf( "%-10s %15s %10s \n", "==========", "MP2 ENERGY", "==========" );
    printf( "%20.16f       \n", total_sum ); 
    printf( "%-10s %15s %10s \n", "==========", "==========", "==========" );
  }
  //
  //
  // SETTING UP SLAVES 
  //
  //
  else{
    while( !done ){
      //cout << "TASK " << taskid << " RECEIVING SOME WORK M8!! " << flush; 
      MPI_Recv( &which_work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
      //cout << "WORK = " << which_work << endl; 
      if( status.MPI_TAG == DIE_TAG /*|| which_work == nwork*/ ){
        //cout << "TASK " << taskid << " JK M8, DYING....!! " << endl; 
        break;
      }
      //
      //
      // creating the actual diagonal part
      //
      //
      psia = which_work + nocc;
      gamma_moint_to_vector( dble_arr1, indx_arr1, corr1matr, evecsXd, norb,
         arr_size1, psia, (psia+1), 0, nocc, 0, nocc, nocc, (psia+1), -10 ); 
      sorter( dble_arr1, indx_arr1, arr_size1 );
      gamma_moint_to_vector( dble_arr2, indx_arr2, corr2matr, evecsXd, norb,
         arr_size2, psia, (psia+1), 0, nocc, 0, nocc, nocc, (psia+1), -10 ); 
      sorter( dble_arr2, indx_arr2, arr_size2 );
      temp_sum = 0.0;
      for( psii = 0; psii < nocc; ++psii  ){
      for( psij = 0; psij < nocc; ++psij ){
      for( psib = nocc; psib <= psia; ++psib ){ 
          denom = evals[ psii ] + evals[ psij ] - evals[ psia ] - evals[ psib ];
          moint1 = get_moint( dble_arr1, indx_arr1, arr_size1, psii, psia, psij, psib );
          moint2 = get_moint( dble_arr2, indx_arr2, arr_size2, psij, psia, psii, psib );
          if( psia == psib )
            temp_sum += ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
          else
            temp_sum += 2. * ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
          //cout << "FOR PSIA, PSII : " << psia << ", " << psii << endl;
          //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psii, psii, psia, moint1 );
          //printf( "(%3d %3d | %3d %3d ) = %20.16f \n", psia, psia, psii, psii, moint2 );
      }
      }
      }
      //cout << "TASK " << taskid << " DONE WITH WORK, SENDING BACK TO 0!! " << endl; 
      MPI_Send( &which_work, 1, MPI_INT, 0, 0, MPI_COMM_WORLD ); 
      MPI_Send( &temp_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD ); 
    } 
  }
}

//    for( int i = 0; i < nocc; ++i ){
//      for( int j = 0; j < nocc; ++j ){
//        for( int b = nocc; b <= a; ++b ){ 
//          moint1 = 0.0;
//          moint2 = 0.0;
//          denom = eigenvals[ i ] + eigenvals[ j ] - eigenvals[ a ] - eigenvals[ b ];
//          moindex1 = get_ijkl( i, a, j, b ); 
//          //
//          //
//          // NEEDS TO BE ORDERED FOR BINARY SEARCH
//          //
//          //
//          if( sort_mo ){
//            mo_index_in_arr = binary_search( indx_arr, moindex1, (size_t)0, size_arr-1 ); 
//            if( mo_index_in_arr >= 0 ){
//              moint1 = dble_arr[ mo_index_in_arr ];
//            }else{
//              //if( taskid == 7 ){ it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex1 );
//              //if( it != ( indx_arr.begin() + size_arr ) ){
//              //  cout << "YA GOOFED" << endl;
//              //  exit( EXIT_FAILURE );
//              //}
//              //}
//            }
//          }else{
//            it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex1 );
//            if( it != ( indx_arr.begin() + size_arr ) ){ 
//              mo_index_in_arr = it - indx_arr.begin(); 
//              moint1 = dble_arr[ mo_index_in_arr ];
//            }
//          }
//          moindex2 = get_ijkl( j, a, i, b ); 
//          //
//          //
//          // NEEDS TO BE ORDERED FOR BINARY SEARCH
//          //
//          //
//          if( sort_mo ){
//            mo_index_in_arr = binary_search( indx_arr, moindex2, (size_t)0, size_arr-1 ); 
//            if( mo_index_in_arr >= 0 ){
//              moint2 = dble_arr[ mo_index_in_arr ];
//            }
//          }else{
//            it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex2 );
//            if( it != ( indx_arr.begin() + size_arr ) ){ 
//              mo_index_in_arr = it - indx_arr.begin(); 
//              moint2 = dble_arr[ mo_index_in_arr ];
//            }
//          }
//          if( a == b )
//            mp2en += ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
//          else
//            mp2en += 2. * ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
//          //if( taskid == 7 ) cout << "TASK 7 : " << moint1 << " " << moint2 << " " << denom << " " << eigenvals[ j ] << " " << \
//                            eigenvals[ i ] << " " << eigenvals[ a ] << " " << eigenvals[ b ] << " " << mp2en << endl;
//        }
//      }
//    }
//    t = clock() - t;
//    double time_other = t / (double)CLOCKS_PER_SEC;
//    //cout << " o OTHER MP2 TIME : " << time_other << endl;
//    //cout << " o TOTAL MP2 TIME : " << ( time_other + time_sort ) << endl;
//  } // end loop over 'a', and so we have this contribution now to mp2 energy
//
//  cout << "TASK " << taskid << " MP2 ENERGY : " << setw( 20 ) << setprecision( 16 ) << mp2en << endl;
//  MPI_Reduce( &mp2en, &total_mp2en, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
