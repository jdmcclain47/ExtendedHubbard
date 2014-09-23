#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <list>
#include <stdio.h>
#include "options.h"
#include "parser.h"
#include "ppp.h"
#include "ewald.h"
#include "cellinfo.h"
#include "ao_ints.h"
#include "Eigen/Dense"


Eigen::Matrix3d make_inv_matr(
    SuperCell& SCell
){
    Eigen::Matrix3d invmat;
    invmat.row(0) = ( SCell.T.row(1) ).cross( (SCell.T.row(2) )) / ( SCell.T.row( 0 ).dot( ( SCell.T.row(1) ).cross( (SCell.T.row(2) ))) );
    invmat.row(1) = ( SCell.T.row(2) ).cross( (SCell.T.row(0) )) / ( SCell.T.row( 0 ).dot( ( SCell.T.row(1) ).cross( (SCell.T.row(2) ))) ); 
    invmat.row(2) = ( SCell.T.row(0) ).cross( (SCell.T.row(1) )) / ( SCell.T.row( 0 ).dot( ( SCell.T.row(1) ).cross( (SCell.T.row(2) ))) );
    return invmat;
}

void inversionSym::makeUniqueList( std::vector< Eigen::Vector3i >& tList, int nx, int ny, int nz ){
  printf( "Looking for inversion symmetry for supercell lattice \n" );
  printf( "   o Periodicity defined after (%3d, %3d, %3d) translations... \n\n", nx, ny, nz );
  int i0, i1, i2, j0, j1, j2;
  for( std::vector< Eigen::Vector3i >::iterator iter_i = tList.begin(); iter_i != tList.end(); ++iter_i ){
    i0 = (*iter_i)( 0 ) % nx;
    i1 = (*iter_i)( 1 ) % ny;
    i2 = (*iter_i)( 2 ) % nz;
    i0 = ( i0 >= 0 ) ? i0
                     : i0 + nx;
    i1 = ( i1 >= 0 ) ? i1
                     : i1 + ny;
    i2 = ( i2 >= 0 ) ? i2
                     : i2 + nz;
    if( invij[ iter_i - tList.begin() ] < 0 ){ // then we haven't found the inverted translational vector
      for( std::vector< Eigen::Vector3i >::iterator iter_j = iter_i; iter_j != tList.end(); ++iter_j ){
        j0 = ( - (*iter_j)( 0 ) ) % nx;
        j1 = ( - (*iter_j)( 1 ) ) % ny;
        j2 = ( - (*iter_j)( 2 ) ) % nz;
        j0 = ( j0 >= 0 ) ? j0
                         : j0 + nx;
        j1 = ( j1 >= 0 ) ? j1
                         : j1 + ny;
        j2 = ( j2 >= 0 ) ? j2
                         : j2 + nz;
        
        if( i0 == j0 &&
            i1 == j1 &&
            i2 == j2 ){
          j0 = (*iter_j)( 0 );
          j1 = (*iter_j)( 1 );
          j2 = (*iter_j)( 2 );
          //printf( "(%3d, %3d, %3d) is inversion of (%3d, %3d, %3d) \n", i0, i1, i2, j0, j1, j2 );
          int indexi = iter_i - tList.begin();
          int indexj = iter_j - tList.begin();
          invij[ indexi ] = indexj; 
          if( indexi != indexj ){ // setting transpose if i!=j
            //printf( "(%3d, %3d, %3d) is inversion of (%3d, %3d, %3d) \n", j0, j1, j2, i0, i1, i2 );
            invij[ indexj ] = indexi; 
          }
          break;
        }
      }
    }
  }
}

void inversionSym::setup( UnitCell& UCell, SuperCell& SCell ){
    naoSuperCell = SCell.nao;
    naoUnitCell = UCell.nao;
    invij.resize( naoSuperCell );
    for( int i = 0; i < naoSuperCell; ++i ) invij[ i ] = -1; // initializing to -1
    makeUniqueList( SCell.reduced_t, SCell.nkx, SCell.nky, SCell.nkz );
}

int aoIntegralFactory::getPerElement( int iao, int jao, int which_cell ){
    int index;
    index = ( iao > jao ) ? ( iao * ( iao + 1 ) ) / 2 + jao + inversion.invij[ which_cell ] * naoUnit_sqp1
                          : ( jao * ( jao + 1 ) ) / 2 + iao + which_cell * naoUnit_sqp1;
    return index; 
}

void aoIntegralFactory::printPerSuperMatrix( 
    SuperCell& SCell, 
    std::vector< double >& inmatr, 
    const char* title
){
    printPerSuperMatrix( SCell, inmatr, title, 12, 8 );
}

void aoIntegralFactory::printPerSuperMatrix( 
    SuperCell& SCell, 
    std::vector< double >& inmatr, 
    const char* title,
    int width,
    int precision
){
    int icount = 0;
    int jcount = 0;
    bool found;
    printf( "------------------------------------------------\n" );
    printf( "%s \n", title );
    printf( "------------------------------------------------\n" );
    for( int itrans = 0; itrans < nTrans; ++itrans, icount++ ){
    for( int i = 0; i < naoUnitCell; ++i ){
      jcount = 0;
      for( int jtrans = 0; jtrans < nTrans; ++jtrans, jcount++ ){
      for( int j = 0; j < naoUnitCell; ++j ){
        Eigen::Vector3i tvec = SCell.reduced_t[ jtrans ] - SCell.reduced_t[ itrans ];
        int index = SCell.get_supercell_index( tvec( 0 ), tvec( 1 ), tvec( 2 ), &found );
        printf( "%*.*f ", width, precision, inmatr[ getPerElement( i, j, index ) ] );
        //printf( "ATOM_I:%3d   R_I:(%3d,%3d,%3d)  ATOM_J:%3d   R_J:(%3d,%3d,%3d)   R_J-R_I(INDEX=%3d):(%3d,%3d,%3d)   DIST:%20.16f \n",
        //        i,SCell.reduced_t[itrans](0),SCell.reduced_t[itrans](1),SCell.reduced_t[itrans](2),
        //        j,SCell.reduced_t[jtrans](0),SCell.reduced_t[jtrans](1),SCell.reduced_t[jtrans](2),
        //        index,SCell.reduced_t[index](0),SCell.reduced_t[index](1),SCell.reduced_t[index](2),
        //        inmatr[ getPerElement( i, j, index ) ] );
      }
      }
      printf( "\n");
    }
    }
}

void aoIntegralFactory::setDefaultOptions(){

    printf( "Setting default options..." );

    Opts.AddOptionString( "XCMETHOD", "WS" );
    Opts.AddOptionString( "COULOMBMETHOD", "EWALD" );
    Opts.AddOptionString( "PPPKERNEL", "COULOMB" );
    Opts.AddOptionDouble( "PPPKERNEL_DIST", 0.0 );
    Opts.AddOptionBool( "PRINT_MATRICES", true );
    Opts.AddOptionBool( "READ_CHECKPOINT", true );
    Opts.AddOptionInt( "EWALD_TOL", -5 );

    /* the following should be copied after we read in the values from the option file */
    SetXCKernel(      Opts.GetOptionString( "XCMETHOD" )      );
    SetCoulombKernel( Opts.GetOptionString( "COULOMBMETHOD" ) );
    SetPPPKernel(     Opts.GetOptionString( "PPPKERNEL" )     );
    pppkerneldist =      Opts.GetOptionDouble( "PPPKERNEL_DIST" );
    tolEwald = Opts.GetOptionInt( "EWALD_TOL" );
    printMatr = Opts.GetOptionBool( "PRINT_MATRICES" );
    readChkpt = Opts.GetOptionBool( "READ_CHECKPOINT" );

    printf( "done!\n" );
}

void aoIntegralFactory::Init( UnitCell& UCell, SuperCell& SCell, const char* ao_options_filename, const char* ppp_options_filename ){
    /* first setting default options */
    setDefaultOptions();
    /* setting the various options */
    std::vector< std::string > options_str = ParseText( ao_options_filename );
    Opts.SetOptionsFromList( options_str );
    /* now setting member values based on these options */
    SetXCKernel(      Opts.GetOptionString( "XCMETHOD" )      );
    SetCoulombKernel( Opts.GetOptionString( "COULOMBMETHOD" ) );
    SetPPPKernel(     Opts.GetOptionString( "PPPKERNEL" )     );
    pppkerneldist =      Opts.GetOptionDouble( "PPPKERNEL_DIST" );
    tolEwald = Opts.GetOptionInt( "EWALD_TOL" );
    printMatr = Opts.GetOptionBool( "PRINT_MATRICES" );
    readChkpt = Opts.GetOptionBool( "READ_CHECKPOINT" );
    /* adding an extra option based on input values */
    do_xc_ppp_correction      = ( pppkerneldist > 0.0 ) && ( pppkern != COULOMB ) && ( xckern      != WS ) && ( xckern      != NONE );
    do_coulomb_ppp_correction = ( pppkerneldist > 0.0 ) && ( pppkern != COULOMB ) && ( coulombkern != WS ) && ( coulombkern != NONE );
    do_ppp_kernel_correction  = (do_xc_ppp_correction || do_coulomb_ppp_correction ); 


    /* setting some needed values */
    inversion.setup( UCell, SCell );
    dim = UCell.dim;
    naoSuperCell = SCell.nao;
    naoUnitCell = UCell.nao;
    naoUnit_sqp1 = ( naoUnitCell * ( naoUnitCell + 1 ) ) / 2;
    nTrans = SCell.total_number_of_cells;
    /* a matrix useful in doing minimum image convention */
    invmat = make_inv_matr( SCell );
    /* making atomic distance matrix */
    makeDistMatr( UCell, SCell, true );
    /* a matrix useful in constructing the full gamma point matrix */
    make_scell_matr( SCell );

    std::string ewaldFile = ".EWALD_DATA";
    bool ewald_read_attempt = false;
    if( readChkpt )
      ewald_read_attempt = read_ewald( UCell, SCell, tolEwald, ewaldFile );
    if( !readChkpt || !ewald_read_attempt ){
      /* now we set up the ewald potential, self potential, and nuclear energy */
      ewald_self = self_potential_ewald_converge_FAST( SCell.T, SCell.K, tolEwald );
      /* by calculating the self potential with this tolerance, the other ewald potentials should be converged as well */
      makeEwaldMatr( UCell, SCell, tolEwald );
      //nucnuc = energy_ewald_converge( UCell.T, UCell.K, UCell.coords, UCell.charge, tolEwald );
      nucnuc = energy_ewald_converge_FAST( UCell.T, UCell.K, UCell.coords, UCell.charge, tolEwald );
    }
    printf( "EWALD SELF ENERGY    : %20.16f \n", ewald_self );
    printf( "EWALD NUCLEAR ENERGY : %20.16f \n", nucnuc );
    write_ewald( UCell, SCell, tolEwald, ewaldFile );

    /* now setting PPP parameters */
    if( strcmp( ao_options_filename, ppp_options_filename ) == 0  ) PPP.Init( options_str, ppp_options_filename );
    else                                                            PPP.Init( ppp_options_filename );
    /* making hopping matrix */
    makeHoppingMatr( UCell, SCell );
    /* making ppp integrals for cut-off */
    makeCutMatrWS( UCell, SCell );
    /* getting prefactor for integrals */
    prefac = PPP.get_prefac();
    /* adding in a correction if we use a periodic 1/r type operator */
    if( do_ppp_kernel_correction ) {
      makeNonCoulombMatr( UCell, SCell, pppkerneldist );
    }
    /* correcting the nuclear-nuclear energy, NOTE : we include terms i = i since we interact with replicated images */
    if( do_ppp_kernel_correction ) {
      for( int i = 0; i < naoUnitCell; ++i ){
      for( int j = i; j < naoUnitCell; ++j ){
        for( int itrans = 0; itrans < nTrans; ++itrans ){
          int index = getPerElement( i, j, itrans );
          nucnuc += aoNonCoulombMatr[ index ];
        }
      }
      }
    }

    /* now we print the options */
    Opts.PrintOptions();
    if( printMatr ){
      printPerSuperMatrix( SCell, aoDistMatr, "ATOMIC DISTANCES", 20, 16 );
      printPerSuperMatrix( SCell, aoHopMatr, "HOPPING MATRIX", 20, 16 );
      printPerSuperMatrix( SCell, aoEwaldMatr, "EWALD MATRIX", 20, 16 );
      printPerSuperMatrix( SCell, aoCutMatrWS, "WS CUT-OFF MATRIX", 20, 16 );
      if( do_ppp_kernel_correction )
        printPerSuperMatrix( SCell, aoNonCoulombMatr, "PPP NON-COULOMB CORRECTION" );
    }
}
