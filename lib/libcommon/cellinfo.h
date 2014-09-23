#ifndef CELLINFO_H
#define CELLINFO_H

#include <string>
#include <math.h>
#include <vector>
#include <cassert>
#include "Eigen/Dense"

static std::vector<std::string> element_names;
static double angs2au = 1.88972598858;

typedef Eigen::Matrix<double,1,3> V3d;


struct UnitCell {

   bool in_angstrom;
   double volume; //volume of unit cell
   int dim;  //dimension of system
   int isym;
   int nat;
   int nao;
   int nShells;
   int nelecs;
   int nOcc;
   int nattype;
   
   std::vector< int > attypelist;

   std::vector< int    > orbitals;
   std::vector< double > charge;
   std::vector< int    > Element;
   std::vector< std::string > ElementStr;
   Eigen::Matrix3d T;  //Translational vectors
   Eigen::Matrix3d K;  //Reciprocal vectors
   std::vector< Eigen::Vector3d > coords;  //vector of coordinates


   void Init(){ Init( "POSCAR" ); };
   void Init( const char* ifilename );
   void ReadXYZFile(const char* xyzFile);
   void pseudo1D_add_vacuum( double vacuum ){
     assert( std::fabs( T( 1, 0 ) ) +
             std::fabs( T( 1, 2 ) ) +
             std::fabs( T( 2, 0 ) ) +
             std::fabs( T( 2, 1 ) ) < 1e-15 
             && "ONLY SET UP FOR THE SIMPLE e1,e2,e3 BASIS VECTORS" );
     T( 1, 1 ) += vacuum; 
     T( 2, 2 ) += vacuum; 
     makeKandVolume();
   };

   private:
      void AngsToAU(){
         T = T * angs2au;
         for( int i = 0; i < nat; ++i )
             coords[ i ] = coords[ i ] * angs2au;
      }
      void WriteOut( );
      void AddAtom(int Element, V3d Atom_Coords);
      void CreatePreDefinedStructure(const char* StructureName);
      void makeKandVolume(); 
};

struct SuperCell {

   /* the total number of unit cells in the supercell
    */
   int total_number_of_cells;
   int nkpt;
   int isym;
   int nao;
   int nkx, nky, nkz;
   double volume;

   int dim;
   /* Defines the translational vectors of the unitcell
    * to generate the supercell
    */
   std::vector< Eigen::Vector3d >  translations;
   std::vector< Eigen::Vector3i >  reduced_t;
   std::vector< Eigen::Vector3d >  kpoints; 
   std::vector< Eigen::Vector3i >  reduced_k; 
   std::vector< Eigen::Vector3d >  coords; 
   std::vector< int    > Element;
   std::vector< std::string > ElementStr;

   /* Gives which unit cell and which atom in the reference 
    * cell the coordinate refers to
    * (unit cell #, atom #)
    */
   Eigen::Matrix<int,Eigen::Dynamic,2> which_ucell;

   /* Defines the translational vectors of supercell
    */
   Eigen::Matrix3d T;
   Eigen::Matrix3d K;

   void Init(UnitCell&);
   void create_supercell_xyz_file( const char* filename );
   int get_supercell_index( int i, int j, int k, bool* found );
   int get_kpoint_index( int i, int j, int k, bool* found );
   Eigen::RowVector3i get_reduced_t_from_index( int index, bool* found );

   private:
      void WriteOutSC();
};

std::string EIntToString(int atomtype);
int EStringToInt(std::string atomtype);
void init_elements();
int gcd(int,int);

#endif
