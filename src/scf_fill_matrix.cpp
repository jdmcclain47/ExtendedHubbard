#include <vector>
#include <iostream>
#include "scf.h"
#include "Eigen/Dense"
#include "cellinfo.h"
#include "common.h"

using namespace std;
using namespace Eigen;

void SCF::Create_Gamma_MatrixXd( 
   UnitCell& UCell, 
   SuperCell& SCell,
   VMatrixXd& block_matrix,
   VMatrixXd& gamma_matrix 
){
   int which_block;
   Vector3i blockij;
   for( int i = 0; i < SCell.total_number_of_cells; i++ ){
      for( int j = i; j < SCell.total_number_of_cells; j++ ){
         if( UCell.dim == 0 ){
            which_block = 0;
         }
         if( UCell.dim == 1 ){
            which_block = j - i;
            if( which_block >= SCell.nkx ){
               which_block -= SCell.nkx;
            }
            if( which_block < 0 ){
               which_block += SCell.nkx;
            }
         }
         if( UCell.dim == 2 ){
            blockij = SCell.reduced_t[ j ] - SCell.reduced_t[ i ]; 
            if( blockij( 0 ) >= SCell.nkx ) blockij( 0 ) -= SCell.nkx;
            if( blockij( 1 ) >= SCell.nky ) blockij( 1 ) -= SCell.nky;

            if( blockij( 0 ) < 0 ) blockij( 0 ) += SCell.nkx;
            if( blockij( 1 ) < 0 ) blockij( 1 ) += SCell.nky;

            which_block = blockij( 0 ) + SCell.nkx * blockij( 1 );
         } 
         if( UCell.dim == 3 ){
            blockij = SCell.reduced_t[ j ] - SCell.reduced_t[ i ]; 
            if( blockij( 0 ) >= SCell.nkx ) blockij( 0 ) -= SCell.nkx;
            if( blockij( 1 ) >= SCell.nky ) blockij( 1 ) -= SCell.nky;
            if( blockij( 2 ) >= SCell.nkz ) blockij( 2 ) -= SCell.nkz;

            if( blockij( 0 ) < 0 ) blockij( 0 ) += SCell.nkx;
            if( blockij( 1 ) < 0 ) blockij( 1 ) += SCell.nky;
            if( blockij( 2 ) < 0 ) blockij( 2 ) += SCell.nkz;

            which_block = blockij( 0 ) + SCell.nkx * blockij( 1 ) + SCell.nkx * SCell.nky * blockij( 2 );
         } 
         gamma_matrix.irrep(0).block( i*nmo_ucell, j*nmo_ucell, nmo_ucell, nmo_ucell ) = block_matrix.irrep( which_block );
      }
      // ... now that we have the upper diagonal blocks of the block matrix, we set the lower triangular portion equal to the upper transpose
      for( int l = 0; l < i; ++l ){
         gamma_matrix.irrep(0).block( i*nmo_ucell, l*nmo_ucell, nmo_ucell, nmo_ucell ) = \
                                gamma_matrix.irrep( 0 ).block( l*nmo_ucell, i*nmo_ucell, nmo_ucell, nmo_ucell ).transpose();
      }
   }

   bool symmetric = true;
   for( int i = 0; i < nmo_scell; ++i ){
       for( int j = 0; j < i; ++j ){
           if( fabs( gamma_matrix.irrep( 0 )( i, j ) - gamma_matrix.irrep( 0 )( j, i ) ) > 1e-08 )
               symmetric = false;
       }
   }
   if( !symmetric ){
      cout << "GAMMA MATRIX IS NOT SYMMETRIC... SOMETHINGS MESSED UP" << endl;
      //exit( EXIT_FAILURE );
   }
}

void SCF::Create_Block_MatrixXd( 
   UnitCell& UCell, 
   SuperCell& SCell,
   VMatrixXd& block_matrix,
   VMatrixXd& gamma_matrix 
){
   for( int ir = 0; ir < nirreps; ++ir ){
     for( int i = 0; i < nmo_ucell; ++i ) 
     for( int j = 0; j < nmo_ucell; ++j )
     block_matrix.irrep( ir )( i, j ) = gamma_matrix.irrep( 0 ).block( 0, ir*nmo_ucell, nmo_ucell, nmo_ucell )( i, j );
   }
}

void SCF::Create_Block_MatrixXd( 
   UnitCell& UCell, 
   SuperCell& SCell,
   VMatrixXd& block_matrix,
   MatrixXd& gamma_matrix 
){
   for( int ir = 0; ir < nirreps; ++ir ){
     for( int i = 0; i < nmo_ucell; ++i ) 
     for( int j = 0; j < nmo_ucell; ++j )
     block_matrix.irrep( ir )( i, j ) = gamma_matrix.block( 0, ir*nmo_ucell, nmo_ucell, nmo_ucell )( i, j );
   }
}
