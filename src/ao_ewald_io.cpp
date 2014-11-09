#include <cstring>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>
#include "cellinfo.h"
#include "ao_ints.h"

using namespace std;

bool aoIntegralFactory::read_ewald( UnitCell& UCell, SuperCell& SCell, int tolEwald, std::string filename ){
    int read_in_nkpt, read_in_nkx, read_in_nky, read_in_nkz;
    int read_in_nao_ucell;
    int read_in_tolEwald;
    vector< double > in_coords;
    vector< double > in_trans;
    bool failed = false;
    FILE* iFile;

    cout << "ATTEMPTING TO READ IN EWALD FROM '" << filename << "'..." << endl;

    iFile = fopen( filename.c_str(), "rb" );
    long int file_size;
    if( iFile == NULL ) failed = true;
    if( failed ){ 
      cout << "FAILED : file doesn't exist." << endl; 
      return false;
    }
    fseek( iFile, 0, SEEK_END );
    file_size = ftell( iFile ); 
    if( file_size == 0 ) failed = true;
    if( failed ){ 
      cout << "FAILED : file empty." << endl; 
      fclose( iFile ); 
      return false;
    }
    rewind( iFile );

    //
    // reading in number of ao's in unit cell 
    //
    fread( &read_in_nao_ucell, sizeof( int ), 1, iFile );
    if( UCell.nao != read_in_nao_ucell ){ failed = true; };
    if( failed ){
      cout << "FAILED : read_in_nao_ucell = " << read_in_nao_ucell << \
              ", current nao = "              << UCell.nao << endl;
      fclose( iFile );
      return false;
    }

    // 
    // reading in coordinates for unit cell
    //
    in_coords.resize( 3 * read_in_nao_ucell );
    fread( &in_coords[ 0 ], sizeof( double ), in_coords.size(), iFile );
    for( int i = 0; i < read_in_nao_ucell; ++i ){
      double coordx = UCell.coords[ i ]( 0 );
      double coordy = UCell.coords[ i ]( 1 );
      double coordz = UCell.coords[ i ]( 2 );
      if( in_coords[ 3 * i + 0 ] != coordx ) failed = true;
      if( in_coords[ 3 * i + 1 ] != coordy ) failed = true;
      if( in_coords[ 3 * i + 2 ] != coordz ) failed = true;
    } 
    if( failed ){
      cout << "FAILED : read_in_coords = " << endl;
      for( int i = 0; i < read_in_nao_ucell; ++i ) printf( "%12.8f %12.8f %12.8f \n",
               in_coords[ 3 * i + 0 ],in_coords[ 3 * i + 1 ],in_coords[ 3 * i + 2 ] );
      cout << "FAILED : current coords = " << endl;
      for( int i = 0; i < read_in_nao_ucell; ++i ) printf( "%12.8f %12.8f %12.8f \n",
               UCell.coords[ i ]( 0 ),UCell.coords[ i ]( 1 ),UCell.coords[ i ]( 2 ) );
      fclose( iFile );
      return false;
    }

    // 
    // reading in translation vectors
    // 
    in_trans.resize( 9 );
    fread( &in_trans[ 0 ], sizeof( double ), in_trans.size(), iFile );
    for( int i = 0; i < 3; ++i ){
      double coordx = UCell.T( i , 0 );
      double coordy = UCell.T( i , 1 );
      double coordz = UCell.T( i , 2 );
      if( in_trans[ 3 * i + 0 ] != coordx ) failed = true;
      if( in_trans[ 3 * i + 1 ] != coordy ) failed = true;
      if( in_trans[ 3 * i + 2 ] != coordz ) failed = true;
    } 
    if( failed ){
      cout << "FAILED : read_in_trans = " << endl;
      for( int i = 0; i < 3; ++i ) printf( "%12.8f %12.8f %12.8f \n",
               in_trans[ 3 * i + 0 ],in_trans[ 3 * i + 1 ],in_trans[ 3 * i + 2 ] );
      cout << "FAILED : current trans = " << endl;
      for( int i = 0; i < 3; ++i ) printf( "%12.8f %12.8f %12.8f \n",
               UCell.T( i , 0 ),UCell.T( i , 1 ),UCell.T( i , 2 ) );
      fclose( iFile );
      return false;
    }

    // 
    // reading in kpoints 
    //
    fread( &read_in_nkpt, sizeof( int ), 1, iFile );
    fread( &read_in_nkx, sizeof( int ), 1, iFile );
    fread( &read_in_nky, sizeof( int ), 1, iFile );
    fread( &read_in_nkz, sizeof( int ), 1, iFile );
    if( read_in_nkpt != SCell.nkpt ) failed = true;
    if( read_in_nkx != SCell.nkx ) failed = true;
    if( read_in_nky != SCell.nky ) failed = true;
    if( read_in_nkz != SCell.nkz ) failed = true;
    if( failed ){
      cout << "FAILED : read_in_nkpt = " << read_in_nkpt << ", current nkpt = " << SCell.nkpt << endl;
      cout << "FAILED : read_in_nkx  = " << read_in_nkx  << ", current nkx  = " << SCell.nkx  << endl;
      cout << "FAILED : read_in_nky  = " << read_in_nky  << ", current nky  = " << SCell.nky  << endl;
      cout << "FAILED : read_in_nkz  = " << read_in_nkz  << ", current nkz  = " << SCell.nkz  << endl;
      fclose( iFile );
      return false;
    }

    // 
    // reading in ewald tolerance 
    //
    fread( &read_in_tolEwald, sizeof( int ), 1, iFile );
    if( read_in_tolEwald != tolEwald ) failed = true;
    if( failed ){
      cout << "FAILED : read_in_tolEwald = " << read_in_tolEwald << ", current tolEwald = " << tolEwald << endl;
      fclose( iFile );
      return false;
    }

    cout << " o SYSTEM CHECK OK! READING EWALD... " << endl;
    fread( &nucnuc, sizeof( double ), 1, iFile );
    fread( &ewald_self, sizeof( double ), 1, iFile );
     
    int nEwald;
    fread( &nEwald, sizeof( int ), 1, iFile );
    if( aoEwaldMatr.empty() )
      aoEwaldMatr.resize( nEwald );
    fread( &aoEwaldMatr[ 0 ], sizeof( double ), nEwald, iFile );
    fclose( iFile );

    return true;
}

void aoIntegralFactory::write_ewald( UnitCell& UCell, SuperCell& SCell, int tolEwald, std::string filename ){
    int read_in_nkpt, read_in_nkx, read_in_nky, read_in_nkz;
    int read_in_nao_ucell;
    int read_in_tolEwald;
    vector< double > in_coords;
    FILE* oFile;

    cout << "WRITING EWALD INFO TO FILE '" << filename << "'..." << endl;

    oFile = fopen( filename.c_str(), "wb" );
    fwrite( &UCell.nao, sizeof( int ), 1, oFile );
    for( int i = 0; i < UCell.nao; ++i ){
      double coordx = UCell.coords[ i ]( 0 );
      double coordy = UCell.coords[ i ]( 1 );
      double coordz = UCell.coords[ i ]( 2 );
      fwrite( &coordx, sizeof( double ), 1, oFile );
      fwrite( &coordy, sizeof( double ), 1, oFile );
      fwrite( &coordz, sizeof( double ), 1, oFile );
    }
    for( int i = 0; i < 3; ++i ){
      double coordx = UCell.T( i , 0 );
      double coordy = UCell.T( i , 1 );
      double coordz = UCell.T( i , 2 );
      fwrite( &coordx, sizeof( double ), 1, oFile );
      fwrite( &coordy, sizeof( double ), 1, oFile );
      fwrite( &coordz, sizeof( double ), 1, oFile );
    }
    fwrite( &SCell.nkpt, sizeof( int ), 1, oFile );
    fwrite( &SCell.nkx,  sizeof( int ), 1, oFile );
    fwrite( &SCell.nky,  sizeof( int ), 1, oFile );
    fwrite( &SCell.nkz,  sizeof( int ), 1, oFile );
    fwrite( &tolEwald,   sizeof( int ), 1, oFile );
    fwrite( &nucnuc,     sizeof( double ), 1, oFile );
    fwrite( &ewald_self, sizeof( double ), 1, oFile );
    int nEwald = aoEwaldMatr.size();
    fwrite( &nEwald,     sizeof( int ), 1, oFile );
    fwrite( &aoEwaldMatr[ 0 ], sizeof( double ), nEwald, oFile );

    fclose( oFile );
}
