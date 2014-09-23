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
#include "moint_read_file.h"
#include "Eigen/Dense"
#include <time.h>
using namespace std;
using namespace Eigen;

void read_gamma_mointb_ind_p( 
    std::string ofilename,
    const int& p_index,
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    size_t& size_arr,
    bool& optimize_size_arr
){
    int i, j, k, l;
    size_t ij_index, kl_index, mo_index;
    double val;
    FILE* oFile;

    bool create_new_kl;
    bool create_first_cont,
         create_second_cont,
         create_third_cont,
         create_fourth_cont;
    std::string filename;
    std::stringstream ss;
    ss << ofilename << p_index;
    ss >> filename;
    oFile = fopen( filename.c_str(), "rb" );

    long int file_size;
    fseek( oFile, 0, SEEK_END );
    file_size = ftell( oFile ); 
    rewind( oFile );

    int header_size = 4 * sizeof( int ) + sizeof( double );
    size_arr = ( file_size / header_size );

    if( optimize_size_arr ){
      //std::cout << "READING FILE '" << filename << "'..." << std::endl;
      //std::cout << "SIZE OF FILE   : " << file_size << std::endl;
      //std::cout << "SIZE OF ENTRY  : " << header_size << std::endl;
      //std::cout << "NELEMENTS      : " << ( file_size / header_size ) << std::endl;
      return;
    }

    for( size_t ind = 0; ind < size_arr; ++ind ){
      fread( &i, sizeof( int ), 1, oFile );
      fread( &j, sizeof( int ), 1, oFile );
      fread( &k, sizeof( int ), 1, oFile );
      fread( &l, sizeof( int ), 1, oFile );
      fread( &val, sizeof( double ), 1, oFile );
      ij_index = ( i * ( i + 1 ) ) / 2 + j;
      kl_index = ( k * ( k + 1 ) ) / 2 + l;
      mo_index = ( ij_index * ( ij_index + 1 ) ) / 2 + kl_index;
      dble_arr[ ind ] = val;
      indx_arr[ ind ] = mo_index;
    } 

    fclose( oFile ); 
}
