#include "Eigen/Dense"
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
#include "common.h"
#include "cellinfo.h"
using namespace std;
using namespace Eigen;

void SCF::write_binary_density( UnitCell& UCell, SuperCell& SCell ){
          write_binary_density( UCell, SCell, ".DENSITY" );
}

void SCF::write_binary_density( UnitCell& UCell, SuperCell& SCell, const char* density_file ){
   int irrep_size = UCell.nao * UCell.nao * sizeof( double );
   int tvec_size = 3 * sizeof( int );
   ofstream outfile( density_file, ios::binary );

   int symmetry = 0;
   outfile.write( (char*)&symmetry, sizeof( symmetry ) );
   outfile.write( (char*)&UCell.nao, sizeof( UCell.nao ) );
   outfile.write( (char*)&nirreps, sizeof( nirreps ) );
   bool found;
   RowVector3i tvec;
   for( int i = 0; i < nirreps; ++i ){
     tvec = SCell.get_reduced_t_from_index( i, &found );
     outfile.write( (char*) tvec.data(), tvec_size );
     outfile.write( (char*) DensXd.irrep( i ).data(), irrep_size ); 
   }

   outfile.close();
}

bool SCF::read_binary_density( UnitCell& UCell, SuperCell& SCell ){
          read_binary_density( UCell, SCell, ".DENSITY" );
}

static bool is_empty( std::ifstream& pFile )
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

bool SCF::read_binary_density( UnitCell& UCell, SuperCell& SCell, const char* density_file ){
   int irrep_size = UCell.nao * UCell.nao * sizeof( double );
   int tvec_size = 3 * sizeof( int );
   ifstream infile( density_file, ios::binary );
   if( !infile.is_open() ){
      printf( "ERROR : read_binary_density, could not read density matrix from file '%s' \n", density_file );
      return false;
   }
   if( is_empty( infile ) ){
      printf( "WARNING : read_binary_density, file '%s' is empty. \n", density_file );
      return false;
   }

   int in_symmetry;
   int in_nmo_ucell;
   int in_nirreps;
   int index;
   printf( "Attempting to read in density from file '%s' \n", density_file ); 
   infile.read( (char*)(&in_symmetry), sizeof( in_symmetry ) );
   infile.read( (char*)(&in_nmo_ucell), sizeof( int ) );
   infile.read( (char*)(&in_nirreps), sizeof( int ) );
   printf( "  - Number of ao's in unit cell : %3d \n", in_nmo_ucell ); 

   bool found;
   RowVector3i tvec;
   if( in_symmetry == 0 && in_nmo_ucell == UCell.nao ){
       MatrixXd temp_dens;
       temp_dens.resize( UCell.nao, UCell.nao );
       for( int i = 0; i < in_nirreps; ++i ){
         infile.read( (char*) tvec.data(), tvec_size );
         index = SCell.get_supercell_index( tvec( 0 ), tvec( 1 ), tvec( 2 ), &found );
         infile.read( (char*) temp_dens.data(), irrep_size );
         if( found ){
           DensXd.irrep( index ) = temp_dens;
         }
       }
   }

   infile.close(); 
   printf( "  - Successfully read! \n" );
   return true;
}
