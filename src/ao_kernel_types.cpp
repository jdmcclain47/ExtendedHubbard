#include "ao_ints.h"
#include <string.h>
#include <cstdlib>
#include <stdio.h>

static double au2angs = 0.529177249;
static double samesite = 1e-3;

void aoIntegralFactory::SetPPPKernel( const char* inkern ){
  if( strcmp( inkern, "MN" ) == 0 ){
    pppkern = MN;
    return;
  }
  if( strcmp( inkern, "OHNO" ) == 0 ){
    pppkern = MN;
    return;
  }
  if( strcmp( inkern, "COULOMB" ) == 0 ){
    pppkern = COULOMB;
    return;
  }
  if( strcmp( inkern, "MAZUMDAR" ) == 0 ){
    pppkern = MAZUMDAR;
    return;
  }
  printf( "ERROR : PPP Kernel '%s' not found \n", inkern );
  exit( EXIT_FAILURE );
}

std::string aoIntegralFactory::getPPPKernelStr(){
  std::string kernmn       ( "MN" );
  std::string kernohno     ( "OHNO" );
  std::string kerncoulomb  ( "COULOMB" );
  std::string kernmazumdar ( "MAZUMDAR" );
  if( pppkern == MN )       return kernmn;
  if( pppkern == OHNO )     return kernohno;
  if( pppkern == COULOMB )  return kerncoulomb; 
  if( pppkern == MAZUMDAR ) return kernmazumdar;
}

double aoIntegralFactory::getPPPInt( 
  std::string ao1, 
  std::string ao2, 
  double dist
){
  return aoIntegralFactory::getPPPInt( ao1, ao2, dist, pppkern );
}

double aoIntegralFactory::getPPPInt( 
  std::string ao1, 
  std::string ao2, 
  double dist,
  const pppkernel inkern
){
  static double mazfac = 0.6117 * au2angs * au2angs;
  if( inkern == MN ){
    double aval = 2. / ( PPP.get_hubbard_u( ao1 ) + PPP.get_hubbard_u( ao2 ) );
    return ( 1./( dist + aval ) );
  }
  if( inkern == OHNO ){
    double aval = 2. / ( PPP.get_hubbard_u( ao1 ) + PPP.get_hubbard_u( ao2 ) );
    return ( 1./sqrt( dist*dist + aval*aval) );
  }
  if( inkern == COULOMB ){
    if( dist < samesite ){ printf( "WARNING : Calling coulomb integral with distance %14.10f... setting to zero. \n", dist ); return 0.0; }
    return ( 1./( dist ) );
  }
  // note : for mazumdar, the 0.6117 is just the 1/14.397 (which is just eV2hartree/au2angs) multiplied
  //        by the usual U of 11.26 eV for conjugated polymers
  //
  // The Ohno potential is usually written in the form :
  //
  //                           U
  //           ---------------------------------
  //           sqrt[ 1 + ( U r_ij / 14.397 )^2 ]
  //
  //
  if( inkern == MAZUMDAR ){
    if( dist > samesite )
      return ( ( ( PPP.get_hubbard_u( ao1 ) + PPP.get_hubbard_u( ao2 ) ) / 2. ) / ( PPP.get_diel() * sqrt( 1. + mazfac * dist * dist ) ) );
    else
      return ( ( ( PPP.get_hubbard_u( ao1 ) + PPP.get_hubbard_u( ao2 ) ) / 2. ) ); // dist = 0 and diel = 1, so we can get rid of the denominator
  }
}


void aoIntegralFactory::SetXCKernel( const char* inkern ){
  setKernel( xckern, inkern );
}

void aoIntegralFactory::SetCoulombKernel( const char* inkern ){
  setKernel( coulombkern, inkern );
}

std::string aoIntegralFactory::getKernelStr( const kernel& inkern ){
  std::string wskern( "WS" );
  std::string ewaldkern( "EWALD" );
  std::string ewaldmsekern( "EWALDMSE" );
  std::string nonekern( "NONE" );
  if( inkern == WS ) return wskern;
  if( inkern == EWALD ) return ewaldkern;
  if( inkern == EWALDMSE ) return ewaldmsekern;
  if( inkern == NONE ) return nonekern;
}

void aoIntegralFactory::setKernel( kernel& inkern, const char* kernName ){
  if( strcmp( kernName, "WS" ) == 0 ){
    inkern = WS;
    return;
  }
  if( strcmp( kernName, "EWALD" ) == 0 ){
    inkern = EWALD; 
    return;
  }
  if( strcmp( kernName, "EWALDMSE" ) == 0 ){
    inkern = EWALDMSE;
    return;
  }
  if( strcmp( kernName, "NONE" ) == 0 ){
    inkern = NONE;
    return;
  }
  printf( "ERROR : Kernel '%s' not found \n", kernName );
  exit( EXIT_FAILURE );
}

double aoIntegralFactory::getXCInt( int iat, int jat, int which_cell ){
  assert( iat < naoUnitCell && jat < naoUnitCell && "iat and jat must be in the unit cell!" );
  int inversion_index = getPerElement( iat, jat, which_cell );
  double val = 0.0;
  if( do_xc_ppp_correction ) val = aoNonCoulombMatr[ inversion_index ];
  if( xckern == WS ) return (prefac*( aoCutMatrWS[ inversion_index ] + val ));
  if( xckern == EWALD ) return (prefac *( aoEwaldMatr[ inversion_index ] + val ));
  if( xckern == EWALDMSE ) return (prefac*( aoEwaldMatr[ inversion_index ] - ewald_self + val ));
  if( xckern == NONE ) return 0.0;
  printf( "ERROR : getXCInt not set up for xckern = %d (see aoints.h for the details of the enumerated kernel type) \n", xckern );
  exit( EXIT_FAILURE );
}


double aoIntegralFactory::getCoulombInt( int iat, int jat, int which_cell ){
  assert( iat < naoUnitCell && jat < naoUnitCell && "iat and jat must be in the unit cell!" );
  int inversion_index = getPerElement( iat, jat, which_cell );
  double val = 0.0;
  if( do_coulomb_ppp_correction ) val = aoNonCoulombMatr[ inversion_index ];
  if( coulombkern == WS ) return (prefac*( aoCutMatrWS[ inversion_index ] + val ) );
  if( coulombkern == EWALD ) return (prefac*( aoEwaldMatr[ inversion_index ] + val ));
  if( coulombkern == EWALDMSE ) return (prefac*( aoEwaldMatr[ inversion_index ] - ewald_self + val ));
  if( coulombkern == NONE ) return 0.0;
  printf( "ERROR : getCoulombInt not set up for xckern = %d (see aoints.h for the details of the enumerated kernel type) \n", xckern );
  exit( EXIT_FAILURE );
}
