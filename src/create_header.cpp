#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include "ao_ints.h"
#include "common.h"
#include "cellinfo.h"
#include "create_header.h"
#include <cstring>

std::string create_kernel_suffix(
  UnitCell& UCell,
  SuperCell& SCell,
  aoIntegralFactory& aoints
){
  // can do UCell and SCell stuff later... but really we just wanna
  // make sure what kernels the integrals and coefficients come from
  std::string header( "_PK_" );
  header.append( aoints.getPPPKernelStr() );
  std::string xc( "_XK_" );
  header.append( xc );
  header.append( aoints.getXCKernelStr() );
  std::string coulomb( "_CK_" ); 
  header.append( coulomb );
  header.append( aoints.getCoulombKernelStr() );
  return header;
}

void write_to_buffer( char* buffer, const size_t max_size, size_t& current_size, std::string& instr ){
  size_t length = instr.length();
  assert( length + 1 + current_size < max_size && "buffer needs a bigger length to write to" ); // +1 to account for \0
  for( int i = 0; i < length; ++i ) buffer[ i + current_size ] = instr[ i ];
  buffer[ length + current_size ] = '\0';
  current_size += (length + 1);
}

char* make_kernel_header(
  UnitCell& UCell,
  SuperCell& SCell,
  aoIntegralFactory& aoints,
  size_t& size
){
  size_t max_size = 100;
  size_t current_size = 0;
  char buffer[ max_size ];

  std::string pppk = aoints.getPPPKernelStr();
  write_to_buffer( buffer, max_size, current_size, pppk );
  std::string xck = aoints.getXCKernelStr();
  write_to_buffer( buffer, max_size, current_size, xck );
  std::string couk = aoints.getCoulombKernelStr(); 
  write_to_buffer( buffer, max_size, current_size, couk );

  size = current_size;

//  std::ofstream output( "TEMP.OUT", std::ios::binary ); 
//  output.write( (char*) &current_size, sizeof( size_t ) );
//  output.write( buffer, current_size );
//  output.close();
  
}

char* read_kernel_header( const char* filename ){
  size_t max_size = 100;
  size_t current_loc = 0;

  std::ifstream input( filename, std::ios::binary );
  input.read( (char*) &max_size, sizeof( size_t ) );
  char* buffer = new char[ max_size ];
  input.read( buffer, max_size );

  while( current_loc < max_size ){
    size_t len = 0;
    while( buffer[ current_loc ] != '\0' ){ len++; current_loc++; }
    if( current_loc < max_size ){
      std::string str( buffer, current_loc - len, len );
    }
    current_loc++; // starting after \0
  }

  delete[] buffer;
}
/*
void check_kernel_header(
  UnitCell& UCell,
  SuperCell& SCell,
  aoIntegralFactory& aoints,
  size_t& size,
  const char* filename 
){
  char* current_header = make_kernel_header( UCell, SCell, aoints, size );
  char* input_header = read_kernel_header( filename );
}
*/

//void header_reader( std::ifstream& infile, const char* format )
