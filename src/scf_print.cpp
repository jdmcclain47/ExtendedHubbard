#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include "scf.h"
#include "Eigen/Dense"
#include "cellinfo.h"
using namespace std;
using namespace Eigen;

void SCF::printPartialMatrix( const char* title, VMatrixXd& in_matr ){
    printf( "Matrix ---- %-15s \n", title );
    int val_per_line = 10;
    int current_col = 0;
    int current_row = 0;
    int offset = 0;
    int vp1 = val_per_line + 1;
    while( offset < nmo_scell ){
      /* first we print "          1            2      ...  " if we're on the first column and first row */
      if( current_row == 0 ){
        for( int i = offset; i < min( nmo_scell+1, offset + vp1 ); ++i ){
          if( i == offset )
            printf( "%-4c", ' ' );
          else
            printf( "%14d", i );
        } 
        printf( "\n" );
        fflush( stdout );
        current_row++;
      }
      else{
        for( int i = offset; i < min( nmo_scell+1, offset + vp1 ); ++i ){
          int which_irrep;
          int which_j; 
          if( i > 0 ){
            which_irrep = (int)floor( (i-1) / nmo_ucell );
            which_j     = (i-1) - which_irrep * nmo_ucell;
          }
          if( i == offset )
            printf( "%-4d", current_row ); 
          else
            printf( "%14.10f", in_matr.irrep( which_irrep )( (current_row-1), (which_j) ) );
          fflush( stdout );
        }
        current_row++;
        printf( "\n" );
        fflush( stdout );
      }
      if( current_row % (nmo_ucell+1) == 0 ){
        current_row -= (nmo_ucell+1); // restarting back to 0
        offset += val_per_line;  
        printf( "\n\n" );
      } 
    }
}
