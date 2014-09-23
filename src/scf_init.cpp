#include "Eigen/Dense"
#include <vector>
#include <cstring>
#include <complex>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "scf.h"
#include "common.h"
#include "parser.h"

using namespace std;
using namespace Eigen;

void SCF::setDefaultOptions(){
    printf( "Setting default scf options... \n" );

    Opts.AddOptionInt( "EN_TOL", -8 );
    Opts.AddOptionInt( "DENS_TOL", -8 );
    Opts.AddOptionInt( "EIGEN_TOL", -6 );
    Opts.AddOptionInt( "NDIIS", 20 );
    Opts.AddOptionInt( "MAX_SCF", 40 );
    Opts.AddOptionBool( "PRINT_MATRICES", true );
    Opts.AddOptionBool( "READ_DENSITY", true );
    Opts.AddOptionBool( "READ_MO", true );
    Opts.AddOptionBool( "USE_DIIS", true );
    Opts.AddOptionBool( "CONVERGE_PSEUDO1D", false );
    Opts.AddOptionBool( "CONVERGE_PSEUDO2D", false );
    Opts.AddOptionString( "PBC", "GAMMA" );

    /* the following should be copied after we read in the values from the option file */
    eTol = pow( 10., Opts.GetOptionInt( "EN_TOL" ) );
    dTol = pow( 10., Opts.GetOptionInt( "DENS_TOL" ) );
    evTol = pow( 10., Opts.GetOptionInt( "EIGEN_TOL" ) );
    ndiis = Opts.GetOptionInt( "NDIIS" );
    max_scf_iter = Opts.GetOptionInt( "MAX_SCF" );
    write_matrices = Opts.GetOptionBool( "PRINT_MATRICES" );
    readin_density = Opts.GetOptionBool( "READ_DENSITY" );
    readin_mo = Opts.GetOptionBool( "READ_MO" );
    use_diis = Opts.GetOptionBool( "USE_DIIS" );
    converge_pseudo_1d = Opts.GetOptionBool( "CONVERGE_PSEUDO1D" );
    converge_pseudo_2d = Opts.GetOptionBool( "CONVERGE_PSEUDO2D" );
    pbc( Opts.GetOptionString( "PBC" ) );

    printf( "  - Complete!\n" );
}

void SCF::Init( UnitCell& UCell, SuperCell& SCell ){
    // ... setting the number of AOs and occupied orbitals based on PPP
    nmo_ucell = UCell.nao; 
    nmo_scell = SCell.nao; 

    nirreps = SCell.nkpt;
    nocc = int( nmo_scell / 2.0 ); 

    /* BEGIN : SCF OPTIONS */
    // setting up default options
    setDefaultOptions();
    // getting user-defined options
    vector< string > scfstr = ParseText( "INCAR" );
    Opts.SetOptionsFromList( scfstr );
    eTol = pow( 10., Opts.GetOptionInt( "EN_TOL" ) );
    dTol = pow( 10., Opts.GetOptionInt( "DENS_TOL" ) );
    evTol = pow( 10., Opts.GetOptionInt( "EIGEN_TOL" ) );
    ndiis = Opts.GetOptionInt( "NDIIS" );
    max_scf_iter = Opts.GetOptionInt( "MAX_SCF" );
    write_matrices = Opts.GetOptionBool( "PRINT_MATRICES" );
    readin_density = Opts.GetOptionBool( "READ_DENSITY" );
    readin_mo = Opts.GetOptionBool( "READ_MO" );
    use_diis = Opts.GetOptionBool( "USE_DIIS" );
    converge_pseudo_1d = Opts.GetOptionBool( "CONVERGE_PSEUDO1D" );
    converge_pseudo_2d = Opts.GetOptionBool( "CONVERGE_PSEUDO2D" );
    pbc( Opts.GetOptionString( "PBC" ) );
    Opts.PrintOptions();
    /* END : SCF OPTIONS */

    /* BEGIN : SETTING UP MATRICES FOR GAMMA AND KPOINT */
    CmatXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    SmatXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    EmatXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    TmatXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    OnsiteXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    FockXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    VenXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    DensXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    ZeroDensXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    old_DensXd.Vresize(nirreps,nmo_ucell,nmo_ucell);
    /* END : SETTING UP MATRICES FOR GAMMA AND KPOINT */

    if( use_diis ) DIIS.Init( ndiis, 1 );
    if( pbc.ret_type() == GAMMA ){
       cout << "Performing gamma calculation..." << endl;
       /* Need to construct a #supercell * #(AOs in unit cell)
        * for the gamma calculation.  This is composed of the 
        * "smaller" matrices above, giving the interaction of the
        * reference cell with various other cells
        */
       SCell_FockXd.Vresize(1,nmo_scell,nmo_scell);
       SCell_DensXd.Vresize(1,nmo_scell,nmo_scell);
       e_vecsXd.Vresize(1,nmo_scell,nmo_scell);
       e_vals.Vresize(1,nmo_scell,1);
    }
    else if( pbc.ret_type() == KPOINT ){
       /* Need to allocate the real portions of the fock and density
        * matrix for fourier transforms between real and k-space
        */
       cout << "Performing k-point calculation..." << endl;
       FockXcd.Vresize(nirreps,nmo_ucell,nmo_ucell);
       DensXcd.Vresize(nirreps,nmo_ucell,nmo_ucell);
       SCell_FockXd.Vresize(1,nmo_scell,nmo_scell);
       SCell_DensXd.Vresize(1,nmo_scell,nmo_scell);
       e_vecsXcd.Vresize(nirreps,nmo_ucell,nmo_ucell);
       e_vals.Vresize(nirreps,nmo_ucell,1);
    }
    else{
       cout << "Invalid Jobtype." << endl;
    }
}
