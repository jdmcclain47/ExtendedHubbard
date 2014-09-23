#include "mpi.h"
#include <cassert>
#include <time.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "mointegrals.h"
#include <math.h>
#include "sorter.h"
#include <string>
#include "Eigen/Dense"
#include "common.h"
#include "moint_read_file.h"
#include "binary_search.h"
#include "mp2.h"
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

double mp2_gamma_full_moints(
    UnitCell& UCell,
    SuperCell& SCell,
    VMatrixXd& evals,
    std::vector< double >& moints
){
    double mp2en = 0.0;
    int nocc = (int)(SCell.nao/2);
    int norb = SCell.nao;
    double denom;
    size_t moindex1, moindex2;
    for( int i = 0; i < nocc; ++i ){ 
    for( int j = 0; j < nocc; ++j ){ 
      for( int a = nocc; a < norb; ++a ){ 
      for( int b = nocc; b < norb; ++b ){ 
        denom = evals.irrep( 0 )( i, 0 ) + evals.irrep( 0 )( j, 0 ) - evals.irrep( 0 )( a, 0 ) - evals.irrep( 0 )( b, 0 );
        moindex1 = get_ijkl( i, a, j, b ); 
        moindex2 = get_ijkl( j, a, i, b ); 
        mp2en += ( moints[ moindex1 ] * ( 2. * moints[ moindex1 ] - moints[ moindex2 ] ) ) / denom;
      }
      }
    }
    }
    return mp2en;
}



double mp2_gamma_ind_p(
  UnitCell& UCell,
  SuperCell& SCell,
  VMatrixXd& evals,
  std::string mointfile,
  bool read_in_from_file,
  moIntegralFactory& moint_class,
  int argc, char *argv[]
){
  //std::cout << "INITIALIZING MPI...." << std::endl;
  int numtasks, taskid, source;
  //MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
  std::cout << "MP2 : RUNNING " << numtasks << " TASKS WITH MPI." << std::endl;
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  std::cout << "MP2 : TASK " << taskid << " STARTED!" << std::endl;
  MPI_Status status;

  int tag1, tag2, tag3, tag4, tag5;
  tag1 = 1;
  tag2 = 2;
  tag3 = 3;
  tag4 = 4;
  tag5 = 5;
  double mp2en = 0.0;
  double total_mp2en = 0.0;
  std::vector< double > eigenvals;
  int nocc, norb;
  Eigen::MatrixXd kernel_matr;
  Eigen::MatrixXd evecsXd;
  std::vector< double > kernel_matr_data;
  std::vector< double > evecsXd_data;
  if( taskid == 0 ){
    nocc = (int)(SCell.nao/2);
    norb = SCell.nao;
    eigenvals.resize( norb );
    for( int i = 0; i < norb; ++i ) eigenvals[ i ] = evals.irrep( 0 )( i, 0 );
    kernel_matr = moint_class.get_kernel_matr();
    evecsXd = moint_class.get_evecsXd();

    /* Send these quantities to all tasks */
    for( int dest = 1; dest < numtasks; ++dest ){
      cout << "SENDING DATA TO TASK " << taskid << endl;
      MPI_Send( &nocc, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD );
      MPI_Send( &norb, 1, MPI_INT, dest, tag2, MPI_COMM_WORLD );
      MPI_Send( &eigenvals[ 0 ], norb, MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD );
      //for( int i = 0; i < norb; ++i ) cout << "TASK " << taskid << " EIGENVALUE " << i << " : " << eigenvals[ i ] << endl;
      MPI_Send( kernel_matr.data(), norb*norb, MPI_DOUBLE, dest, tag4, MPI_COMM_WORLD );
      MPI_Send( evecsXd.data(), norb*norb, MPI_DOUBLE, dest, tag5, MPI_COMM_WORLD );
      cout << "DONE : SENDING DATA TO TASK " << taskid << endl;
    }
    //cout << taskid << " KERNEL MATRIX : " << endl;
    //cout << kernel_matr.block( 0, 0, 3, 3 ) << endl;
    //cout << taskid << " EVECS MATRIX  : " << endl;
    //cout << evecsXd.block( 0, 0, 3, 3 ) << endl;
  }
  if( taskid > 0 ){
    for( int dest = 1; dest < numtasks; ++dest ){
      if( taskid == dest ){
        cout << "TASK " << dest << " RECIEVING DATA" << endl;
        MPI_Recv( &nocc, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status );
        MPI_Recv( &norb, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD, &status );
        eigenvals.resize( norb );
        kernel_matr_data.resize( norb * norb );
        evecsXd_data.resize( norb * norb );
        MPI_Recv( &eigenvals[ 0 ], norb, MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD, &status );
        //for( int i = 0; i < norb; ++i ) cout << "TASK " << taskid << " EIGENVALUE " << i << " : " << eigenvals[ i ] << endl;
        MPI_Recv( &kernel_matr_data[ 0 ], norb*norb, MPI_DOUBLE, 0, tag4, MPI_COMM_WORLD, &status );
        //for( int i = 0; i < norb*norb; ++i ) cout << "TASK " << taskid << " KERNEL VAL " << i << " : " << kernel_matr_data[ i ] << endl;
        MPI_Recv( &evecsXd_data[ 0 ], norb*norb, MPI_DOUBLE, 0, tag5, MPI_COMM_WORLD, &status );
        //for( int i = 0; i < norb*norb; ++i ) cout << "TASK " << taskid << " EVEVCS_DATA " << i << " : " << evecsXd_data[ i ] << endl;
        kernel_matr.resize( norb, norb );
        evecsXd.resize( norb, norb );
        int counter = 0;
        for( int i = 0; i < norb; ++i )
        for( int j = 0; j < norb; ++j, counter++ ){
          kernel_matr( j, i ) = kernel_matr_data[ counter ];
          evecsXd( j, i ) = evecsXd_data[ counter ];
        }
        //cout << taskid << " KERNEL MATRIX : " << endl;
        //cout << kernel_matr.block( 0, 0, 3, 3 ) << endl;
        //cout << taskid << " EVECS MATRIX  : " << endl;
        //cout << evecsXd.block( 0, 0, 3, 3 ) << endl;
        cout << "DONE : TASK " << dest << " RECIEVING DATA" << endl;
      }
    }
  }

  MPI_Barrier( MPI_COMM_WORLD );

  int mo_index_in_arr;
  double denom, moint1, moint2;
  size_t moindex1, moindex2;
  std::vector< size_t >::iterator it;
  bool sort_mo = true;

  clock_t t, timer;
  size_t opt_size_arr, size_arr;
  std::vector< size_t > indx_arr;
  std::vector< double > dble_arr;
  bool optimize_size_arr;


  if( read_in_from_file ){
    opt_size_arr = 0;
    optimize_size_arr = true;
    std::cout << "READING FROM MULTIPLE MOINT FILES" << std::endl;
    for( int i = 0; i < SCell.nao; ++i ){
      read_gamma_mointb_ind_p( mointfile, i, dble_arr, indx_arr, size_arr, optimize_size_arr );
      if( size_arr > opt_size_arr ) opt_size_arr = size_arr;
    }
    std::cout << "OPTIMAL SIZE OF ARRAYS " << opt_size_arr << std::endl;
    optimize_size_arr = false;
    dble_arr.resize( opt_size_arr );
    indx_arr.resize( opt_size_arr );
  }else{
    dble_arr.resize( 5000000 );
    indx_arr.resize( 5000000 );
  }

  clock_t total_time;
  if( taskid == 0 )
  total_time = clock();

  std::vector< int > break_up_work;
  break_up_work.resize( numtasks + 1 );
  break_up_work[ 0 ] = nocc;
  break_up_work[ numtasks ] = norb;
  int break_up_load = (int)( ( norb - nocc ) / numtasks );

  for( int i = 1; i < break_up_work.size() - 1; ++i ) break_up_work[ i ] = break_up_work[ i - 1 ] + break_up_load;

  //for( int a = nocc; a < norb; ++a ){ 
  for( int i = 0; i < numtasks; ++i ){
  if( taskid == i ){

  int low = break_up_work[ i ];
  int high = break_up_work[ (i + 1) ];

  for( int a = low; a < high; ++a ){ 
    cout << "MP2 COMPLETED ( " << setw( 3 ) << ( a - nocc + 1 ) << " / " << setw( 3 ) << ( norb - nocc ) << " ) " << endl;
    if( read_in_from_file ){
      read_gamma_mointb_ind_p( mointfile, a, dble_arr, indx_arr, size_arr, optimize_size_arr );
    }else{
      t = clock();
      //moint_class.gamma_moint_to_vector( dble_arr, indx_arr, size_arr, 
      //                                   a, (a+1), 0, nocc, 
      //                                   nocc, (a+1), 0, nocc, -10 );
      gamma_moint_to_vector( dble_arr, indx_arr, kernel_matr, evecsXd, norb, size_arr,
                             a, (a+1), 0, nocc, 
                             nocc, (a+1), 0, nocc, -10 );
      t = clock() - t;
      if( taskid == 7 ){
        cout << " MOINT GENERATION FOR TASK " << taskid << endl;
        cout << " - time        = " << t / (double)CLOCKS_PER_SEC << " seconds." << endl;
        cout << " - # non-zero  = " << size_arr << endl;
        cout << " - 'a' virt val= " << a << endl;
      }
    }
    
    t = clock();
    if( sort_mo ){ sorter( dble_arr, indx_arr, size_arr ); }
    t = clock() - t;
    double time_sort = t / (double)CLOCKS_PER_SEC;
    //if( sort_mo )  cout << " o TIME FOR MOINT SORTING : " << time_sort << endl;
    
    t = clock();
    for( int i = 0; i < nocc; ++i ){
      for( int j = 0; j < nocc; ++j ){
        for( int b = nocc; b <= a; ++b ){ 
          moint1 = 0.0;
          moint2 = 0.0;
          denom = eigenvals[ i ] + eigenvals[ j ] - eigenvals[ a ] - eigenvals[ b ];
          moindex1 = get_ijkl( i, a, j, b ); 
          //
          //
          // NEEDS TO BE ORDERED FOR BINARY SEARCH
          //
          //
          if( sort_mo ){
            mo_index_in_arr = binary_search( indx_arr, moindex1, (size_t)0, size_arr-1 ); 
            if( mo_index_in_arr >= 0 ){
              moint1 = dble_arr[ mo_index_in_arr ];
            }else{
              //if( taskid == 7 ){ it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex1 );
              //if( it != ( indx_arr.begin() + size_arr ) ){
              //  cout << "YA GOOFED" << endl;
              //  exit( EXIT_FAILURE );
              //}
              //}
            }
          }else{
            it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex1 );
            if( it != ( indx_arr.begin() + size_arr ) ){ 
              mo_index_in_arr = it - indx_arr.begin(); 
              moint1 = dble_arr[ mo_index_in_arr ];
            }
          }
          moindex2 = get_ijkl( j, a, i, b ); 
          //
          //
          // NEEDS TO BE ORDERED FOR BINARY SEARCH
          //
          //
          if( sort_mo ){
            mo_index_in_arr = binary_search( indx_arr, moindex2, (size_t)0, size_arr-1 ); 
            if( mo_index_in_arr >= 0 ){
              moint2 = dble_arr[ mo_index_in_arr ];
            }
          }else{
            it = find( indx_arr.begin(), indx_arr.begin() + size_arr, moindex2 );
            if( it != ( indx_arr.begin() + size_arr ) ){ 
              mo_index_in_arr = it - indx_arr.begin(); 
              moint2 = dble_arr[ mo_index_in_arr ];
            }
          }
          if( a == b )
            mp2en += ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
          else
            mp2en += 2. * ( moint1 * ( 2. * moint1 - moint2 ) ) / denom;
          //if( taskid == 7 ) cout << "TASK 7 : " << moint1 << " " << moint2 << " " << denom << " " << eigenvals[ j ] << " " << \
                            eigenvals[ i ] << " " << eigenvals[ a ] << " " << eigenvals[ b ] << " " << mp2en << endl;
        }
      }
    }
    t = clock() - t;
    double time_other = t / (double)CLOCKS_PER_SEC;
    //cout << " o OTHER MP2 TIME : " << time_other << endl;
    //cout << " o TOTAL MP2 TIME : " << ( time_other + time_sort ) << endl;
  } // end loop over 'a', and so we have this contribution now to mp2 energy

  cout << "TASK " << taskid << " MP2 ENERGY : " << setw( 20 ) << setprecision( 16 ) << mp2en << endl;
  MPI_Reduce( &mp2en, &total_mp2en, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  
  } // end if == taskid 
  } // end loop over tasks
  //MPI_Finalize();

  MPI_Barrier( MPI_COMM_WORLD );

  return total_mp2en;
}
