#include "iostream"
#include "cassert"
#include "stdio.h"
#include "stdlib.h"
#include "math.h" 
#include "sstream"
using namespace std;

int main( int argc, char* argv[] ){
  double dTol = 1e-15;
  double sval, lval, old_sval, old_lval, volume;
  double length, alpha, tdt, gdg, area, self;
  int nmax;
  stringstream str;
  assert( argc >= 3 && "NEEDS AT LEAST TWO TRAILING ARGUMENTS" );
  str << argv[ 1 ];
  str >> length;
  str.clear();
  str << argv[ 2 ];
  str >> alpha;

  area = length * length;
  volume = area * length;
  printf( "LENGTH OF BOX : %20.16f \n", length );
  printf( "AREA OF BOX : %20.16f \n", area );
  printf( "VOLUME OF BOX : %20.16f \n", volume );
  printf( "ALPHA VALUE   : %20.16f \n", alpha );

  //
  //
  // Converging short ranged term
  //
  //
  old_sval = 999.9; 
  sval = 0.0;
  nmax = 1;
  printf( "%40s \n", " :: CONVERGING SHORT RANGED :: " );
  printf( "%20s %20s %10s \n", "SVAL", "DIFF", "NSUB" );
  while( fabs( old_sval - sval ) > dTol ){ 
  old_sval = sval;
  sval = 0.0;
  for( int ix = -nmax; ix <= nmax; ++ix ){ 
  for( int iy = -nmax; iy <= nmax; ++iy ){ 

    tdt = length * sqrt( (ix*ix+iy*iy) );
    if( abs( ix ) + abs( iy ) != 0 )
    sval += 0.5 * erfc( alpha * tdt ) / tdt; 
    
  }
  }

  printf( "%20.16f %20.16f %10d \n", sval, fabs( sval - old_sval ), nmax );
  nmax++;
  }
  nmax--;
  nmax--;


  //
  //
  // Converging long ranged term
  //
  //
  old_lval = 999.9; 
  lval = 0.0;
  nmax = 1;
  printf( "%40s \n", " :: CONVERGING LONG RANGED :: " );
  printf( "%20s %20s %10s \n", "LVAL", "DIFF", "NSUB" );
  while( fabs( old_lval - lval ) > dTol ){ 
  old_lval = lval;
  lval = 0.0;
  for( int ix = -nmax; ix <= nmax; ++ix ){ 
  for( int iy = -nmax; iy <= nmax; ++iy ){ 

    gdg = 2. * M_PI / length * sqrt( (double)(ix*ix+iy*iy) );
    if( abs( ix ) + abs( iy ) == 0 )
    lval += 0.5 / area * ( - 2. * M_PI / alpha / sqrt( M_PI ) );
    else
    lval += 0.5 / area * ( M_PI / gdg * ( 2. * erfc( gdg/2./alpha) ) );
    
  }
  }

  printf( "%20.16f %20.16f %10d \n", lval, fabs( lval - old_lval ), nmax );
  nmax++;
  }
  nmax--;
  nmax--;

  self = alpha / sqrt( M_PI );
  printf( "2D TOT.     : %20.16f \n", (lval+sval-self) );
  printf( "2D TOT. * L : %20.16f \n", (lval+sval-self)*length );


  //
  //
  // Converging short ranged term
  //
  //
  old_sval = 999.9; 
  sval = 0.0;
  nmax = 1;
  printf( "%40s \n", " :: CONVERGING SHORT RANGED :: " );
  printf( "%20s %20s %10s \n", "SVAL", "DIFF", "NSUB" );
  while( fabs( old_sval - sval ) > dTol ){ 
  old_sval = sval;
  sval = 0.0;
  for( int ix = -nmax; ix <= nmax; ++ix ){ 
  for( int iy = -nmax; iy <= nmax; ++iy ){ 
  for( int iz = -nmax; iz <= nmax; ++iz ){ 

    tdt = length * sqrt( (ix*ix+iy*iy+iz*iz) );
    if( abs( ix ) + abs( iy ) != 0 )
    sval += 0.5 * erfc( alpha * tdt ) / tdt; 
    
  }
  }
  }

  printf( "%20.16f %20.16f %10d \n", sval, fabs( sval - old_sval ), nmax );
  nmax++;
  }
  nmax--;
  nmax--;


  //
  //
  // Converging long ranged term
  //
  //
  old_lval = 999.9; 
  lval = 0.0;
  nmax = 1;
  printf( "%40s \n", " :: CONVERGING LONG RANGED :: " );
  printf( "%20s %20s %10s \n", "LVAL", "DIFF", "NSUB" );
  while( fabs( old_lval - lval ) > dTol ){ 
  old_lval = lval;
  lval = 0.0;
  for( int ix = -nmax; ix <= nmax; ++ix ){ 
  for( int iy = -nmax; iy <= nmax; ++iy ){ 
  for( int iz = -nmax; iz <= nmax; ++iz ){ 

    gdg = 2. * M_PI / length * sqrt( (double)(ix*ix+iy*iy+iz*iz) );
    if( abs( ix ) + abs( iy ) + abs( iz ) != 0 )
    lval += 0.5 * (4.0*M_PI) / volume * exp( - gdg*gdg/4./alpha/alpha ) / (gdg*gdg);
    
  }
  }
  }

  printf( "%20.16f %20.16f %10d \n", lval, fabs( lval - old_lval ), nmax );
  nmax++;
  }
  nmax--;
  nmax--;

  self = alpha / sqrt( M_PI ) + 0.5 * M_PI / alpha / alpha / volume;
  printf( "3D TOT.     : %20.16f \n", (lval+sval-self) );
  printf( "3D TOT. * L : %20.16f \n", (lval+sval-self)*length );

}
