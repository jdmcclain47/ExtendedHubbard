#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "cellinfo.h"
#include "common.h"
#include "mointegrals.h"
#include "Eigen/Dense"
#include <time.h>
using namespace std;
using namespace Eigen;


void moIntegralFactory::write_gamma_moint_to_fcidump( 
    const char* ofilename,
    int tol
 ){
    int ij_index, kl_index, mo_index;
    double val;
    FILE* oFile;
    printf( "WRITING TO FCIDUMP FILE.... \n" );
    oFile = fopen( ofilename, "w" );

    bool create_new_kl;
    bool create_first_cont,
         create_second_cont,
         create_third_cont,
         create_fourth_cont;
    double dtol = pow( 10., 1. * tol );
    //
    //
    // printing header
    //
    //
    int orbsym = 1;
    int ms2 = 0;
    int isym = 1;
    fprintf( oFile, " &FCI NORB=%3d,NELEC=%3d,MS2=%2d,\n", nmo_scell, nmo_scell, ms2 );
    fprintf( oFile, "  ORBSYM=" );
    for( int i = 0; i < nmo_scell; ++i ) 
      fprintf( oFile, "%1d,", orbsym );
    fprintf( oFile, "\n" );
    fprintf( oFile, "  ISYM=%1d \n", isym );
    fprintf( oFile, " &END\n" );
    
    for( int k = 0; k < nmo_scell; ++k ){
    create_first_cont = true;
    printf( "Performing gamma transform ; Completed (%3d /%3d ) ... \n", (k+1), nmo_scell ); cout << flush;
    for( int l = 0; l <= k; ++l ){
      create_second_cont = true;
      kl_index = ( k * ( k + 1 ) ) / 2 + l;
      create_new_kl = true;
      for( int i = 0; i < nmo_scell; ++i ){
      create_third_cont = true;
      for( int j = 0; j <= i; ++j ){
        create_fourth_cont = true;
        ij_index = ( i * ( i + 1 ) ) / 2 + j;
        if( ij_index >= kl_index ){
          mo_index = ( ij_index * ( ij_index + 1 ) ) / 2 + kl_index;
          val = quarter_gamma_transform( k, l, i, j, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
          if( fabs( val ) > dtol )
            fprintf( oFile, "%24.16e %3d %3d %3d %3d \n", val, (i+1), (j+1), (k+1), (l+1) );

          create_first_cont  = false;
          create_second_cont = false;
          create_third_cont  = false;
          create_fourth_cont = false;
        }
      }
      }
    }
    }
    for( int i = 0; i < nmo_scell; ++i ){
      val = evals.irrep( 0 )( i, 0 );
      fprintf( oFile, "%24.16e %3d %3d %3d %3d \n", val, (i+1), (i+1), 0, 0 );
    }
    val = 0.0;
    fprintf( oFile, "%24.16e %3d %3d %3d %3d \n", val, 0, 0, 0, 0 );
    fprintf( oFile, "%24.16e %3d %3d %3d %3d \n", val, 0, 0, 0, 0 );
    fclose( oFile );
}


void moIntegralFactory::write_gamma_mointb_ind_p( 
    std::string ofilename,
    int tol
){
    int ij_index, kl_index, mo_index;
    double val;
    FILE* oFile;
    printf( "WRITING MOINTEGRAL BINARY FILE.... \n" );

    bool create_new_kl;
    bool create_first_cont,
         create_second_cont,
         create_third_cont,
         create_fourth_cont;
    double dtol = pow( 10., 1. * tol );
    for( int i = 0; i < nmo_scell; ++i ){
      std::string filename;
      std::stringstream ss;
      ss << ofilename << i;
      ss >> filename;
      oFile = fopen( filename.c_str(), "wb" );
      create_first_cont = true;
      printf( "WRITING MOINTEGRAL BINARY FILE '%s'.... ( %4d / %4d ) \n", filename.c_str(), i, nmo_scell ); cout << flush;
      for( int j = 0; j <= i; ++j ){
        create_second_cont = true;
        ij_index = ( i * ( i + 1 ) ) / 2 + j;
        create_new_kl = true;
        for( int k = 0; k < nmo_scell; ++k ){
        create_third_cont = true;
        for( int l = 0; l <= k; ++l ){
          create_fourth_cont = true;
          kl_index = ( k * ( k + 1 ) ) / 2 + l;
//          if( ij_index >= kl_index ){
          if( i >= k ){
            mo_index = ( ij_index * ( ij_index + 1 ) ) / 2 + kl_index;
            val = quarter_gamma_transform( i, j, k, l, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
            if( fabs( val ) > dtol ){
              fwrite( &i, sizeof( int ), 1, oFile );
              fwrite( &j, sizeof( int ), 1, oFile );
              fwrite( &k, sizeof( int ), 1, oFile );
              fwrite( &l, sizeof( int ), 1, oFile );
              fwrite( &val, sizeof( double ), 1, oFile );
            }
     
            create_first_cont  = false;
            create_second_cont = false;
            create_third_cont  = false;
            create_fourth_cont = false;
          }
        }
        }
      }
      fclose( oFile );
    }
}


