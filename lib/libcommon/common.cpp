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
#include "common.h"

using namespace std;
using namespace Eigen;


string trim( string str ){
   for( int i = 0; i < str.length(); i++ ){
      if( str[i] == ' ' ){
         str = str.erase( i, 1 ); 
      }
   }
   return str;
}

void split( string line, const char* delimiter, vector< string >* parsed_string){
   string lincpy = line;
   string delim = delimiter;
   string temp;
   int end_position;
   while( lincpy.find( delimiter ) < lincpy.length() ){
     temp = lincpy.substr( 0, lincpy.find( delimiter ) );
     temp = trim( temp );
     (*parsed_string).push_back( temp );
     end_position = lincpy.find( delimiter ) + delim.length();
     lincpy = lincpy.erase( 0, end_position  );
   }
   temp = trim( lincpy );
   if( temp.length() > 0 )
      (*parsed_string).push_back( temp );
}

static int const ElementSize = 87;
static char const* ElementNames[ElementSize] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
	"P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
	"Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
	"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
	"Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "X" };


string getElementName( int AtomicNumber ){
   if( AtomicNumber == 0 ){
      return "X";
   }
   assert( (AtomicNumber >= 1) && (AtomicNumber <= ElementSize) );
   return ElementNames[AtomicNumber-1];
}



int getElementNumber( string Name )
{
   int ret = -1;
   for(int i=0;i<ElementSize;i++ ){
      if ( strcasecmp( Name.c_str(), ElementNames[i] ) == 0 ){
         ret = i + 1;
         return ret;
      }
   }
   if(strcasecmp(Name.c_str(),"X") == 0 ){
      ret = 0;
      return ret;
   }
   return ret;
}


static int const ElementOccupations = 40;
static int const OccupationNumbers[ElementOccupations] = 
	{ 1, 0,
	  2, 0,
	  2, 1, 0,
	  2, 2, 0,
	  2, 2, 0,  1, 0,
	  2, 2, 0,  2, 0,
	  2, 2, 0,  3, 0,
	  2, 2, 0,  4, 0,
	  2, 2, 0,  5, 0,
	  2, 2, 0,  6, 0 }; 



vector<int> getOccupation( int AtomicNumber ){
   assert( AtomicNumber > 0 && AtomicNumber < 11 && "Occupations currently only set from atomic number 1 to 11" );
   vector<int> vec;
   switch( AtomicNumber ){
      case  1: vec.assign ( &OccupationNumbers[0],  &OccupationNumbers[1] );
      case  2: vec.assign ( &OccupationNumbers[2],  &OccupationNumbers[3] );
      case  3: vec.assign ( &OccupationNumbers[4],  &OccupationNumbers[6] );
      case  4: vec.assign ( &OccupationNumbers[7],  &OccupationNumbers[9] );
      case  5: vec.assign ( &OccupationNumbers[10], &OccupationNumbers[14] );
      case  6: vec.assign ( &OccupationNumbers[15], &OccupationNumbers[19] );
      case  7: vec.assign ( &OccupationNumbers[20], &OccupationNumbers[24] );
      case  8: vec.assign ( &OccupationNumbers[25], &OccupationNumbers[29] );
      case  9: vec.assign ( &OccupationNumbers[30], &OccupationNumbers[34] );
      case 10: vec.assign ( &OccupationNumbers[35], &OccupationNumbers[39] );
   }
   return vec;
}

void VMatrixXd::Vresize(int number_irrep, int rows, int cols){
   vec_matrix.resize(number_irrep);
   for(int i=0;i<number_irrep;i++){
      vec_matrix[i].resize(rows,cols);
      vec_matrix[i] = MatrixXd::Zero(rows,cols);
   }
   m_rows = rows;
   m_cols = cols;
   n_irrep = number_irrep;
}

MatrixXd& VMatrixXd::irrep(int ir_rep){
   MatrixXd& block_ir = vec_matrix[ir_rep];
   return block_ir;
}

void VMatrixXcd::Vresize(int number_irrep, int rows, int cols){
   vec_matrix.resize(number_irrep);
   for(int i=0;i<number_irrep;i++){
      vec_matrix[i].resize(rows,cols);
      vec_matrix[i] = MatrixXcd::Zero(rows,cols);
   }
   m_rows = rows;
   m_cols = cols;
   n_irrep = number_irrep;
}

MatrixXcd& VMatrixXcd::irrep(int ir_rep){
   MatrixXcd& block_ir = vec_matrix[ir_rep];
   return block_ir;
}
