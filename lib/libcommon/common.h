#ifndef COMMON_H
#define COMMON_H

#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

enum pbc_type { GAMMA, KPOINT };

class PBC_CLASS {
    public:
      pbc_type type;
      void operator()( pbc_type in_type ){ type = in_type; };
      void operator()( const char* intype ){
        if( strcmp( intype, "GAMMA" ) == 0 ) type = GAMMA; 
        else if( strcmp( intype, "KPOINT" ) == 0 ) type = KPOINT; 
        else { printf( "ERROR : ( class pbc ) No type '%s' found \n", intype ); exit( EXIT_FAILURE ); }
      };
      std::string str(){
        std::string gamma_str ( "GAMMA" );
        if( type == GAMMA ) return gamma_str; 
        std::string kpoint_str ( "KPOINT" );
        if( type == KPOINT ) return kpoint_str; 
        printf( "ERROR : ( class pbc : str() ) type not initialized? \n" ); exit( EXIT_FAILURE );
      };
      pbc_type ret_type(){ return type; };
};


std::string trim( std::string str );
void split( std::string line, const char* delimiter, std::vector< std::string >* parsed_string);

std::string getElementName( int atomic_number );
int getElementNumber( std::string element_name );

class VMatrixXd {

   int m_rows;
   int m_cols;
   int n_irrep;
   std::vector<Eigen::MatrixXd> vec_matrix;	//a vector of matrices in matrix form
   
   public:
      void Vresize(int number_irrep, int rows, int cols);
      Eigen::MatrixXd& irrep(int);
      int size(){return vec_matrix.size();};
};

class VMatrixXcd {

   int m_rows;
   int m_cols;
   int n_irrep;
   std::vector<Eigen::MatrixXcd> vec_matrix;	//a vector of matrices in matrix form
   
   public:
      void Vresize(int number_irrep, int rows, int cols);
      Eigen::MatrixXcd& irrep(int);
      int size(){return vec_matrix.size();};
};

#endif
