#include "Eigen/Dense"
#include <vector>
#include <cstring>
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

void SCF::write_binary_mo_coeff( 
   UnitCell& UCell, 
   SuperCell& SCell, 
   aoIntegralFactory& aoints
){
   write_binary_mo_coeff( UCell, SCell, aoints, ".MOCOEFF" );
}

void SCF::write_binary_mo_coeff( 
   UnitCell& UCell, 
   SuperCell& SCell, 
   aoIntegralFactory& aoints,
   const char* mocoeff_file 
){
   int irrep_sizeXd = SCell.nao * SCell.nao * sizeof( double );
   int irrep_sizeXcd = UCell.nao * UCell.nao * sizeof( complex< double > );
   int eval_gamma  = SCell.nao * sizeof( double );
   int eval_kpoint = UCell.nao * sizeof( double );
   char kernelstr[ 10 ];
   std::string pppkern     = aoints.getPPPKernelStr();
   std::string xckern      = aoints.getXCKernelStr();
   std::string coulombkern = aoints.getCoulombKernelStr();
   size_t strlen;
   int tvec_size = 3 * sizeof( int );
   ofstream outfile( mocoeff_file, ios::binary );

   printf( "PRINTING MO COEFF TO FILE '%s'.\n", mocoeff_file );
   int symmetry = 0;
   outfile.write( (char*)&symmetry, sizeof( symmetry ) );
   outfile.write( (char*)&UCell.nao, sizeof( UCell.nao ) );
   outfile.write( (char*)&nirreps, sizeof( nirreps ) );
   outfile.write( (char*)&pbc.type, sizeof( pbc.type ) );
   // writing the length of string and the string for kernels
   strcpy( kernelstr, pppkern.c_str() );
   strlen = pppkern.length();
   outfile.write( (char*)&strlen, sizeof( size_t ) );
   outfile.write( (char*)&kernelstr, sizeof( kernelstr ) );
   strcpy( kernelstr, xckern.c_str() );
   strlen = xckern.length();
   outfile.write( (char*)&strlen, sizeof( size_t ) );
   outfile.write( (char*)&kernelstr, sizeof( kernelstr ) );
   strcpy( kernelstr, coulombkern.c_str() );
   strlen = coulombkern.length();
   outfile.write( (char*)&strlen, sizeof( size_t ) );
   outfile.write( (char*)&kernelstr, sizeof( kernelstr ) );
   
   bool found;
   RowVector3i tvec;
   if( pbc.type == GAMMA ){
     outfile.write( (char*) e_vecsXd.irrep( 0 ).data(), irrep_sizeXd ); 
     outfile.write( (char*) e_vals.irrep( 0 ).data(), eval_gamma ); 
   }
   if( pbc.type == KPOINT ){
     for( int i = 0; i < nirreps; ++i ){
       outfile.write( (char*) e_vecsXcd.irrep( i ).data(), irrep_sizeXcd ); 
       outfile.write( (char*) e_vals.irrep( i ).data(), eval_kpoint ); 
     }
   }
    

   outfile.close();
   printf( "  - Complete! \n" );
}

bool SCF::read_binary_mo_coeff(
   UnitCell& UCell, 
   SuperCell& SCell, 
   aoIntegralFactory& aoints
){
   return read_binary_mo_coeff( UCell, SCell, aoints, ".MOCOEFF" );
}

static bool is_empty( std::ifstream& pFile )
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

bool SCF::read_binary_mo_coeff( 
   UnitCell& UCell, 
   SuperCell& SCell, 
   aoIntegralFactory& aoints,
   const char* mocoeff_file 
){
   bool read_success = false;
   int irrep_sizeXd = SCell.nao * SCell.nao * sizeof( double );
   int irrep_sizeXcd = UCell.nao * UCell.nao * sizeof( complex< double > );
   int eval_gamma  = SCell.nao * sizeof( double );
   int eval_kpoint = UCell.nao * sizeof( double );
   char kernelstr[ 10 ];
   int tvec_size = 3 * sizeof( int );
   ifstream infile( mocoeff_file, ios::binary );
   if( !infile.is_open() ) printf( "ERROR : read_binary_mo_coeff, could not read from file '%s' \n", mocoeff_file );
   if( is_empty( infile ) ) printf( "WARNING : read_binary_mo_coeff, file '%s' is empty. \n", mocoeff_file );

   int in_symmetry;
   int in_nmo_ucell;
   int in_nirreps;
   PBC_CLASS in_pbc;
   int index;
   size_t strsize;
   
   printf( "Attempting to read MO coefficients from file '%s' \n", mocoeff_file ); 
   infile.read( (char*)(&in_symmetry), sizeof( in_symmetry ) );
   infile.read( (char*)(&in_nmo_ucell), sizeof( int ) );
   infile.read( (char*)(&in_nirreps), sizeof( int ) );
   infile.read( (char*)(&in_pbc.type), sizeof( in_pbc.type ) );
   // writing the length of string and the string for kernels
   infile.read( (char*)&strsize, sizeof( size_t ) );
   infile.read( (char*)&kernelstr, sizeof( kernelstr ) );
   std::string pppkern( kernelstr, strsize );
   infile.read( (char*)&strsize, sizeof( size_t ) );
   infile.read( (char*)&kernelstr, sizeof( kernelstr ) );
   std::string xckern( kernelstr, strsize );
   infile.read( (char*)&strsize, sizeof( size_t ) );
   infile.read( (char*)&kernelstr, sizeof( kernelstr ) );
   std::string coulombkern( kernelstr, strsize );
          
   printf( "PARAM. FROM FILE          / CURRENT PARAM. \n", mocoeff_file );
   printf( "  - NAO    : %10d   /    %10d \n", in_nmo_ucell, nmo_ucell ); 
   printf( "  - NIRREP : %10d   /    %10d \n", in_nirreps,   nirreps ); 
   printf( "  - PBC    : %10s   /    %10s \n", in_pbc.str().c_str(), pbc.str().c_str() ); 
   printf( "  - PPPK   : %10s   /    %10s \n", pppkern.c_str(), aoints.getPPPKernelStr().c_str() ); 
   printf( "  - XCK    : %10s   /    %10s \n", xckern.c_str(), aoints.getXCKernelStr().c_str() ); 
   printf( "  - COUK   : %10s   /    %10s \n", coulombkern.c_str(), aoints.getCoulombKernelStr().c_str() ); 

   if( in_nmo_ucell == nmo_ucell && 
       in_nirreps   == nirreps   &&
       in_pbc.type  == pbc.type  &&
       pppkern      == aoints.getPPPKernelStr() &&
        xckern      == aoints.getXCKernelStr()  &&
       coulombkern  == aoints.getCoulombKernelStr() 
     ){ 

     bool found;
     RowVector3i tvec;
     if( in_symmetry == 0 && in_nmo_ucell == UCell.nao ){
         MatrixXd temp_matrXd;
         MatrixXcd temp_matrXcd;
         if( in_pbc.type == GAMMA ){
           temp_matrXd.resize( SCell.nao, SCell.nao );
           infile.read( (char*) temp_matrXd.data(), irrep_sizeXd );
           e_vecsXd.irrep( 0 ) = temp_matrXd;
           infile.read( (char*) e_vals.irrep( 0 ).data(), eval_gamma ); 
         }
         if( in_pbc.type == KPOINT ){
           for( int i = 0; i < nirreps; ++i ){
             infile.read( (char*) e_vecsXcd.irrep( i ).data(), irrep_sizeXcd ); 
             infile.read( (char*) e_vals.irrep( i ).data(), eval_kpoint ); 
           }
         }
     }
    
     read_success = true;
   }
   if( read_success ) printf( " -----Successfully read.\n" );
   if( !read_success ) printf( " -----Failed to read.\n" );

   infile.close(); 
   return read_success;
}
