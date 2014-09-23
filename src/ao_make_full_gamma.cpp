#include "cellinfo.h"
#include "ao_ints.h"
#include "Eigen/Dense"
using namespace Eigen;

void aoIntegralFactory::make_scell_matr(
    SuperCell& SCell
){
   int which_cell;
   Vector3i cellij;
   scell_matr.resize( SCell.total_number_of_cells, SCell.total_number_of_cells );

   for( int i = 0; i < SCell.total_number_of_cells; i++ ){
      for( int j = 0; j < SCell.total_number_of_cells; j++ ){
         if( SCell.dim == 0 ){
            which_cell = 0;
         }
         if( SCell.dim == 1 ){
            which_cell = j - i;
            if( which_cell >= SCell.nkx ){
               which_cell -= SCell.nkx;
            }
            if( which_cell < 0 ){
               which_cell += SCell.nkx;
            }
         }
         if( SCell.dim == 2 ){
            cellij = SCell.reduced_t[ j ] - SCell.reduced_t[ i ];
            if( cellij( 0 ) >= SCell.nkx ) cellij( 0 ) -= SCell.nkx;
            if( cellij( 1 ) >= SCell.nky ) cellij( 1 ) -= SCell.nky;

            if( cellij( 0 ) < 0 ) cellij( 0 ) += SCell.nkx;
            if( cellij( 1 ) < 0 ) cellij( 1 ) += SCell.nky;

            which_cell = cellij( 0 ) + SCell.nkx * cellij( 1 );
         }
         if( SCell.dim == 3 ){
            cellij = SCell.reduced_t[ j ] - SCell.reduced_t[ i ];
            if( cellij( 0 ) >= SCell.nkx ) cellij( 0 ) -= SCell.nkx;
            if( cellij( 1 ) >= SCell.nky ) cellij( 1 ) -= SCell.nky;
            if( cellij( 2 ) >= SCell.nkz ) cellij( 2 ) -= SCell.nkz;

            if( cellij( 0 ) < 0 ) cellij( 0 ) += SCell.nkx;
            if( cellij( 1 ) < 0 ) cellij( 1 ) += SCell.nky;
            if( cellij( 2 ) < 0 ) cellij( 2 ) += SCell.nkz;

            which_cell = cellij( 0 ) + SCell.nkx * cellij( 1 ) + SCell.nkx * SCell.nky * cellij( 2 );
         }
         scell_matr( i, j ) = which_cell;
      }
   }
}
