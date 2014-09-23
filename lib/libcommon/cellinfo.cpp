#include <iostream>
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <fstream>
#include <ctype.h>
#include "cellinfo.h"
#include "parser.h"
#include "common.h"
using namespace std;
using namespace Eigen;



void UnitCell::ReadXYZFile( const char* XYZFile ){
    cout << "Creating unit cell... " << endl;
    vector< string > str = ParseText( XYZFile );
    // ... the file should be set up as so :
    //     ( optional 'angstrom' identifier to identify the coordinates are given in angstroms )
    //     ( optional 'symmetry : X' argument, currently set to 1 for a CNT.... )
    //     DIMENSION OF SYSTEM
    //     T1 (if dim > 0 )
    //     T2 (if dim > 1 ) 
    //     T3 (if dim > 2 ) 
    //     ATOM1   XCOORD1   YCOORD1   ZCOORD1    
    //     ATOM1   XCOORD1   YCOORD1   ZCOORD1    
    //     ATOM2   XCOORD2   YCOORD2   ZCOORD2    
    //     
    // ATOM1 is a string representing an atom
    // 
    // ATOM1 is used if there are any pre-defined hopping elements
    // between the atoms, and using 'X' for ATOM1 will make sure that no
    // pre-defined hamiltonian matrix elements will be used between the 
    // various atomic sites.
    const char* angs = "angstrom";
    const char* sym = "symmetry";
    int offset = 0;
    string strtemp; 
    int itemp;
    {
      istringstream iss( str[ 0 ] );
      iss >> strtemp; 
      if( strcmp( strtemp.c_str(), angs ) == 0 ){
          in_angstrom = true;
          offset += 1;
      }
      iss.str(" ");
      iss.clear();

      iss.str( str[ offset ] );
      iss >> strtemp;
      if( strcmp( strtemp.c_str(), sym ) == 0 ){
          offset += 1;
          iss.str( str[ offset ] );
          iss >> itemp;
          if( itemp == 0 ){
            isym = 0;
            cout << "Symmetry set to default.." << endl;
            offset += 1;
          }
          else if( itemp == 1 ){
            isym = 1;
            cout << "Symmetry set to cnt symmetry.." << endl;
            offset += 1;
          }
          else{
            cout << "ERROR : Set symmetry = '" << itemp << "'.  Only symmetry is either 0 (normal periodic calculation..) or 1 (CNT skew symmetry)" << endl;
            exit( EXIT_FAILURE );
          }
      }
      iss.str(" ");
      iss.clear();
     
      iss.str( str[ offset ] ); // ... dimension of system 
      iss >> dim;
      offset += 1;
    }
    double dtemp;
    for( int i = 0; i < dim; i++ ){
      {
        istringstream iss( str[ i * 3 + 0 + offset ] ); // ... translation vector T(i,0) 
        iss >> T(i,0);
      }
      {
        istringstream iss( str[ i * 3 + 1 + offset ] ); // ... translation vector T(i,1) 
        iss >> T(i,1);
      }
      {
        istringstream iss( str[ i * 3 + 2 + offset ] ); // ... translation vector T(i,2) 
        iss >> T(i,2);
      }
    }
    offset += dim * 3;
    int remaining_strings = str.size() - offset;
    int linesize = 4; // ... the number of strings on a line in the XYZ file (atom, x, y, z )
    nat = ( remaining_strings / linesize );
    if( remaining_strings % linesize != 0 ){
        cout << "ERROR : Error reading XYZ file '" << XYZFile << "'. ";
        cout << "Should follow the format : " << endl;
        cout << "ATOM1   XCOORD1   YCOORD1   ZCOORD1" << endl;
        cout << "ATOM2   XCOORD2   YCOORD2   ZCOORD2" << endl;
        cout << "Parsed XYZ file : " << endl;
        for( int i = 0; i < str.size(); ++i ){
            cout << str[ i ] << endl;
        }
        exit( EXIT_FAILURE );
    }
    string c;
    double temp;
    coords.resize( nat );
    charge.resize( nat );
    Element.resize( nat );
    ElementStr.resize( nat );

    bool found_new_atom;

    for( int i = 0; i < nat; ++i ){
      c = str[ i * linesize + offset];
      if( isalpha( c[ 0 ] ) ){
          charge[ i ] = 1.0;
          Element[ i ] = getElementNumber( c );
          ElementStr[ i ] = c; 
          if( Element[ i ] == -1 ){
              cout << "WARNING : NO ELEMENT NAME FOR '" << c << "'." << endl;
              //exit( EXIT_FAILURE );
          }
          found_new_atom = true;
          if( i == 0 ){ // ... we push_back the element into our unique atom_list
              attypelist.push_back( Element[ i ] );
          }else{ // ... we check whether the current element is different from the previous
              for( int iat = 0; iat < attypelist.size(); ++iat ){
                  if( Element[ i ] == attypelist[ iat ] ){ // ... this is equivalent to a previous atom 
                      found_new_atom = false;
                      break;
                  } 
              }
              if( found_new_atom )
                  attypelist.push_back( Element[ i ] );
          } 
      }else{
          cout << "ReadXYZFile only set-up currently to read atoms as characters" << endl;
          cout << "Parsed XYZ file : " << endl;
          for( int i = 0; i < str.size(); ++i ){
              cout << str[ i ] << endl;
          }
          exit( EXIT_FAILURE );
      }
      {
        istringstream iss( str[ i * linesize + 1 + offset ] ); // ... x coordinate value 
        iss >> temp;
        coords[ i ]( 0 ) = temp;
      }
      {
        istringstream iss( str[ i * linesize + 2 + offset ] ); // ... y coordinate value 
        iss >> temp;
        coords[ i ]( 1 ) = temp;
      }
      {
        istringstream iss( str[ i * linesize + 3 + offset ] ); // ... z coordinate value 
        iss >> temp;
        coords[ i ]( 2 ) = temp;
      }
      //cout << ", coordinate : " << coords( i, 0 ) << ", " << coords( i, 1 ) << ", " << coords( i, 2 ) << endl;
    }  
    
    nattype = attypelist.size();

}

/* These are the only things that depend on T */
void UnitCell::makeKandVolume(){
   volume = T.row(0).dot( T.row(1).cross(T.row(2)) );
   K.row(0) = 2.0*M_PI*T.row(1).cross(T.row(2))/volume;
   K.row(1) = 2.0*M_PI*T.row(2).cross(T.row(0))/volume;
   K.row(2) = 2.0*M_PI*T.row(0).cross(T.row(1))/volume;
   volume = fabs( volume );
}

void UnitCell::Init( const char* ifilename ){
   in_angstrom = false;
   // Seeing whether we used a predefined system and writing out system to the xyzfile 
   //if( JOB.usepredefinedsystem == 1 ){
   //   CreatePreDefinedStructure( JOB.predefinedsystem.c_str() );
   //}
   // Reading in the atoms and their coordinates
   ReadXYZFile( ifilename );
   // Converting everything to atomic units
   if( in_angstrom )
       AngsToAU();

   // ... if our dimension is less than 3, it is easiest to create the reciprocal vectors by 
   //     completing the space of R^3, then generating a reciprocal vector space based on that
   if( dim < 3 && dim > 0 ){
      MatrixXd span_space;
      span_space.resize( 3, 3 );
      span_space = MatrixXd::Random(3,3);
      for( int l = 0; l < dim; l++ ){
         span_space.col(l) = T.row(l).transpose();
      } 
      for( int l = dim; l < 3; l++ ){
         V3d temp_vec;
         temp_vec.transpose() = span_space.col(l);
         for(int i=0;i<l;i++){
            temp_vec.transpose() -= span_space.col(i) * \
		( span_space.col(i).dot( span_space.col(l) ) )/( span_space.col(i).dot( span_space.col(i) ) );
         }
         span_space.col(l) = temp_vec.transpose();
      }
      T = span_space.transpose();
   }

   makeKandVolume();

   nao = nat;
   nOcc = (int)( nao / 2 );

   WriteOut();
}


/*
void UnitCell::CreatePreDefinedStructure(const char* StructureName){
   if(strncmp (StructureName, "ZGNR",4)==0){
      sscanf(StructureName,"ZGNR(%d)",&PDStructure.pds_gnr_length);
      PDStructure.System = "ZGNR";
      PDStructure.make_zgnr();
   }
   if(strncmp (StructureName, "AGNR",4)==0){
      sscanf(StructureName,"AGNR(%d)",&PDStructure.pds_gnr_length);
      PDStructure.System = "AGNR";
      PDStructure.make_agnr();
   }
   if(strncmp (StructureName, "CNT",3)==0){
      sscanf(StructureName,"CNT(%d,%d)",&PDStructure.pds_cnt_n,&PDStructure.pds_cnt_m);
      PDStructure.System = "CNT";
      PDStructure.make_cnt();
   }
   if(strncmp (StructureName, "TPA",3)==0){
      PDStructure.System = "TPA";
      PDStructure.make_tpa();
   }
}
*/

void SuperCell::Init( UnitCell & UCell){
   vector< string > str = ParseText( "KPOINT" );
   dim = UCell.dim;
   total_number_of_cells = 0;
   int nx, ny, nz;
   nx = 0; ny = 0; nz = 0;
   total_number_of_cells = 1;

   isym = UCell.isym;

   volume = UCell.volume;
   T.row( 0 ) = UCell.T.row( 0 );
   T.row( 1 ) = UCell.T.row( 1 );
   T.row( 2 ) = UCell.T.row( 2 );

   if( dim > 0 ){
     istringstream iss( str[ 0 ] ); 
     int itemp;
     iss >> itemp;
     total_number_of_cells = itemp;
     nx = itemp;
     T.row(0) = (UCell.T.row(0)) * nx;
     volume *= nx;
     if( nx <= 0 ){
         cout << "If dim > 0, you must have at least 1 K-Point in x-dir."; 
         cout << "Check file KPOINT and dimension of system." << endl;
         exit( EXIT_FAILURE );
     }
   }
   if( dim > 1 ){
     istringstream iss( str[ 1 ] ); 
     int itemp;
     iss >> itemp;
     total_number_of_cells *= itemp;
     ny = itemp;
     T.row(1) = (UCell.T.row(1)) * ny;
     volume *= ny;
     if( ny <= 0 ){
         cout << "If dim > 1, you must have at least 1 K-Point in y-dir."; 
         cout << "Check file KPOINT and dimension of system." << endl;
         exit( EXIT_FAILURE );
     }
   }
   if( dim > 2 ){
     istringstream iss( str[ 2 ] ); 
     int itemp;
     iss >> itemp;
     total_number_of_cells *= itemp;
     nz = itemp;
     T.row(2) = (UCell.T.row(2)) * nz;
     volume *= nz;
     if( nz <= 0 ){
         cout << "If dim > 1, you must have at least 1 K-Point in z-dir."; 
         cout << "Check file KPOINT and dimension of system." << endl;
         exit( EXIT_FAILURE );
     }
   }

   nao = total_number_of_cells * UCell.nao;

   K.row(0) = 2.0*M_PI*T.row(1).cross(T.row(2))/volume;
   K.row(1) = 2.0*M_PI*T.row(2).cross(T.row(0))/volume;
   K.row(2) = 2.0*M_PI*T.row(0).cross(T.row(1))/volume;

   translations.resize( total_number_of_cells );
   reduced_t.resize( total_number_of_cells );
   Vector3i itempv;
   if( dim == 0 ){
      translations[ 0 ]( 0 ) = 0.0;
      translations[ 0 ]( 1 ) = 0.0;
      translations[ 0 ]( 2 ) = 0.0;
      itempv << 0, 0, 0;
      reduced_t[ 0 ] = itempv;
   }
   if( dim == 1 ){
     for(int i = 0; i < nx; i++){
        translations[ i ] =  1.*i*(UCell.T).row(0);
        itempv << i, 0, 0; 
        reduced_t[ i ] = itempv;
     }
   }
   if( dim == 2 ){
     for(int i = 0; i < nx; i++)
     for(int j = 0; j < ny; j++){
        translations[ i + j * nx ] =  1.*i*(UCell.T).row(0) + 1.*j*(UCell.T).row(1);
        itempv << i, j, 0;
        reduced_t[ i + j * nx ] = itempv;
     }
   }
   if( dim == 3 ){
      for(int i = 0; i < nx; i++)
      for(int j = 0; j < ny; j++)
      for(int k = 0; k < nz; k++){
         translations[ i + j * nx + k * nx * ny ] =  1.*i*(UCell.T).row(0) + 1.*j*(UCell.T).row(1) + 1.*k*(UCell.T).row(2);
         itempv << i, j, k;
         reduced_t[ i + j * nx + k * nx * ny ] = itempv;
      }
   }

   coords.resize(UCell.nat * total_number_of_cells);
   Element.resize( coords.size() );
   ElementStr.resize( coords.size() );
   which_ucell.resize(UCell.nat * total_number_of_cells,2);
   for(int j=0;j<total_number_of_cells;j++){
      for(int i=0;i<UCell.nat;i++){
         coords[ j * UCell.nat + i ] = UCell.coords[ i ] + translations[ j ];
         Element[ j * UCell.nat + i ] = UCell.Element[ i ];
         ElementStr[ j * UCell.nat + i ] = UCell.ElementStr[ i ];
         which_ucell(j * UCell.nat + i,0) = j; // ... this atom in the supercell is in the 'j'-th supercell
         which_ucell(j * UCell.nat + i,1) = i; // ... this atom in the supercell is an image of the 'i'-th atom in the original cell
      }
   }


   nkx = nx; nky = ny; nkz = nz;
   kpoints.resize( total_number_of_cells );
   reduced_k.resize( total_number_of_cells );
   nkpt = total_number_of_cells;
   if( dim == 1 ){
      for(int i = 0; i < nx; i++)
        kpoints[ i ] =  (1./nx)*i*(UCell.K).row(0);
   }
   if( dim == 2 ){
     for(int i = 0; i < nx; i++)
     for(int j = 0; j < ny; j++)
        kpoints[ i + j * nx ] =  (1./nx)*i*(UCell.K).row(0) + (1./ny)*j*(UCell.K).row(1);
   }
   if( dim == 3 ){
     for(int i = 0; i < nx; i++)
     for(int j = 0; j < ny; j++)
     for(int k = 0; k < nz; k++){
        kpoints[ i + j * nx + k * nx * ny ] =  (1./nx)*i*(UCell.K).row(0) + (1./ny)*j*(UCell.K).row(1) + (1./nz)*k*(UCell.K).row(2);
        reduced_k[ i + j * nx + k * nx * ny ]( 0 ) = i;
        reduced_k[ i + j * nx + k * nx * ny ]( 1 ) = j;
        reduced_k[ i + j * nx + k * nx * ny ]( 2 ) = k;
     }
   }

   WriteOutSC();
}

int SuperCell::get_kpoint_index( int i, int j, int k, bool* found ){
    i = ( i >= 0 ) ? i : i + nkx;
    j = ( j >= 0 ) ? j : j + nky;
    k = ( k >= 0 ) ? k : k + nkz;
    if( i >= nkx || j >= nky || k >= nkz ){
        (*found) = false;
        return -1;
    }
    (*found) = true;
    return ( i + j * nkx + k * nkx * nky );
}

int SuperCell::get_supercell_index( int i, int j, int k, bool* found ){
    i = ( i >= 0 ) ? i : i + nkx;
    j = ( j >= 0 ) ? j : j + nky;
    k = ( k >= 0 ) ? k : k + nkz;
    if( i >= nkx || j >= nky || k >= nkz ){
        (*found) = false;
        return -1;
    }
    (*found) = true;
    return ( i + j * nkx + k * nkx * nky );
}

RowVector3i SuperCell::get_reduced_t_from_index( int index, bool* found ){
    RowVector3i outvec;
    if( index > total_number_of_cells ){
        (*found) = false;
        outvec( 0 ) = 999; 
        outvec( 1 ) = 999; 
        outvec( 2 ) = 999; 
        return outvec;
    }
    (*found) = true;
    outvec( 0 ) = index % nkx; 
    outvec( 1 ) = ( ( index - outvec( 0 ) ) / nkx ) % nky;
    outvec( 2 ) = ( ( index - outvec( 0 ) - outvec( 1 ) * nkx ) / nkx / nky ); 
    //cout << "CHECKING..." << index << " = INDEX, " << outvec << " = REDUCED_T : " << get_supercell_index( outvec( 0 ), outvec( 1 ), outvec( 2 ), found ) << " = TRUE INDEX" << endl;
    return outvec;
}

void SuperCell::create_supercell_xyz_file( const char* filename ){
    ofstream output( filename );
    for( int i = 0; i < coords.size(); ++i ){
        output << getElementName( SuperCell::Element[ i ] ) << " " << coords[ i ]( 0 ) << " " << coords[ i ]( 1 ) << " " << coords[ i ]( 2 ) << endl;
    }
    output.close();
}


void SuperCell::WriteOutSC(){
   cout << "****************\n" << endl;
   cout << "SUPER CELL INFO \n" << endl;
   for(int i=0;i<dim;i++){
      printf("Super Cell Translational vector [%d] : %20.16f %20.16f %20.16f \n",i,T(i,0),T(i,1),T(i,2));
   }
//   for(int i=0;i<total_number_of_cells;i++){
//      printf("Image Unit Cell / Reduced Coords [%d] : %12.8f %12.8f %12.8f %3d %3d %3d \n",i,translations[i](0),translations[i](1),translations[i](2), reduced_t[ i ]( 0 ), reduced_t[ i ]( 1 ), reduced_t[ i ]( 2 ));   
//   }
   for(int i=0;i<total_number_of_cells;i++){
      printf("K-vectors [%d] : %12.8f %12.8f %12.8f \n",i,kpoints[i](0),kpoints[i](1),kpoints[i](2));   
   }
//   cout << "ATOMIC POSITIONS : " << endl;
//   for( int i = 0; i < coords.size(); i++ ){
//      printf( "%s  %20.16f  %20.16f  %20.16f \n", getElementName( SuperCell::Element[ i ] ).c_str(), coords[ i ]( 0 ), coords[ i ]( 1 ), coords[ i ]( 2 ) );
//   }
   cout << "****************\n" << endl;
   printf("\n");
}


void UnitCell::WriteOut( ){
   cout << "***************\n" << endl;
   cout << "UNIT CELL INFO \n" << endl;
   printf("Number of unique atoms : %3d \n", nattype ); 
   printf("Unique Atom list : ");
   for( int i = 0; i < (nattype-1); ++i ){
       printf("%s, ", getElementName( attypelist[ i ] ).c_str() );
   }
   printf("%s \n", getElementName( attypelist[ nattype-1 ] ).c_str() );

   if( in_angstrom ){
      printf("Coordinates input in angstroms...\n");
   }
   printf("Dimension of system : %d \n",dim);   
   for(int i=0;i<dim;i++){
      printf("T%d : %12.8f %12.8f %12.8f \n",(i+1),T(i,0),T(i,1),T(i,2));   
   }
   for(int i=0;i<dim;i++){
      printf("G%d : %12.8f %12.8f %12.8f \n",(i+1),K(i,0),K(i,1),K(i,2));   
   }
   printf("Coordinates: \n");
   for(int i=0;i<nat;i++){
      printf("   %s   %12.8f %12.8f %12.8f \n",getElementName(Element[i]).c_str(),coords[i](0),coords[i](1),coords[i](2));
   }
   cout << "***************\n" << endl;
   printf("\n");
}

