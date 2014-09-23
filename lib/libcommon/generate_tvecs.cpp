#include <vector>
#include <cstdio>
#include <complex>
#include <math.h>
#include <iostream>
#include <iomanip>
#include "assert.h"
#include "generate_tvecs.h"
#include "Eigen/Dense" 
using namespace std;
using namespace Eigen;

void generate_tvecs(
  int& nvals,
  std::vector< double >& tmatrix,
  int in_xold, int in_xnew,
  int in_yold, int in_ynew,
  int in_zold, int in_znew,
  double in_xoffset,
  double in_yoffset,
  double in_zoffset,
  std::vector< double >& xyzvals,
  std::vector< double >& xyz_sq_vals
){

assert( tmatrix.size() == 9 );

double xoffset, yoffset, zoffset;
int xold, yold, zold;
int xnew, ynew, znew;
xnew = in_xnew;
ynew = in_ynew;
znew = in_znew;

xold = in_xold;
yold = in_yold;
zold = in_zold;

xoffset = in_xoffset;
yoffset = in_yoffset;
zoffset = in_zoffset;

nvals = ( 2 * xnew + 1 ) * ( 2 * ynew + 1 ) * ( 2 * znew + 1 ) - \
        ( 2 * xold + 1 ) * ( 2 * yold + 1 ) * ( 2 * zold + 1 );
if( xold == 0 && yold == 0 && zold == 0 ){ nvals += 1; }
xyzvals.resize( 3 * nvals );
xyz_sq_vals.resize( 3 * nvals );

std::vector< double > t1 ( 3, 0.0 );
std::vector< double > t2 ( 3, 0.0 );
std::vector< double > t3 ( 3, 0.0 );
t1[ 0 ] = tmatrix[ 0 ]; t1[ 1 ] = tmatrix[ 1 ]; t1[ 2 ] = tmatrix[ 2 ];
t2[ 0 ] = tmatrix[ 3 ]; t2[ 1 ] = tmatrix[ 4 ]; t2[ 2 ] = tmatrix[ 5 ];
t3[ 0 ] = tmatrix[ 6 ]; t3[ 1 ] = tmatrix[ 7 ]; t3[ 2 ] = tmatrix[ 8 ];

//printf( " ===== TVECS ===== \n" );
//printf( " T1 - %12.8f  %12.8f  %12.8f \n", t1[ 0 ], t1[ 1 ], t1[ 2 ] );
//printf( " T2 - %12.8f  %12.8f  %12.8f \n", t2[ 0 ], t2[ 1 ], t2[ 2 ] );
//printf( " T3 - %12.8f  %12.8f  %12.8f \n", t3[ 0 ], t3[ 1 ], t3[ 2 ] );

double xv, yv, zv;
double xv2, yv2, zv2;
double x1, x2, x3;
double y1, y2, y3;
double z1, z2, z3;
//
//
// now doing the (0,0,0) translation
//
//
xv = 0.0; xv2 = xv * xv;
yv = 0.0; yv2 = yv * yv;
zv = 0.0; zv2 = zv * zv;
int counter = 0;
if( xold == 0 && yold == 0 && zold == 0 ){
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
}

//
//
// now doing the (x,0,0) type translations, where x > 0
//
//
for( int ix = xold + 1; ix <= xnew; ++ix ){
for( int iy = 0; iy <= 0; ++iy ){
for( int iz = 0; iz <= 0; ++iz ){
  xv = xoffset;
  yv = yoffset;
  zv = zoffset;
  // x > 0 contribution
  xv += 1. * ix * t1[ 0 ];
  yv += 1. * ix * t1[ 1 ];
  zv += 1. * ix * t1[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  // x < 0 contribution
  xv -= 2. * ix * t1[ 0 ];
  yv -= 2. * ix * t1[ 1 ];
  zv -= 2. * ix * t1[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
}
}
}
//
//
// now doing the (0,y,0) type translations, where y > 0
//
//
for( int ix = 0; ix <= 0; ++ix ){
for( int iy = yold + 1; iy <= ynew; ++iy ){
for( int iz = 0; iz <= 0; ++iz ){
  xv = xoffset;
  yv = yoffset;
  zv = zoffset;
  // y > 0 contribution
  xv += 1. * iy * t2[ 0 ];
  yv += 1. * iy * t2[ 1 ];
  zv += 1. * iy * t2[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  // y < 0 contribution
  xv -= 2. * iy * t2[ 0 ];
  yv -= 2. * iy * t2[ 1 ];
  zv -= 2. * iy * t2[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
}
}
}
for( int ix = 0; ix <= 0; ++ix ){
for( int iy = 0; iy <= 0; ++iy ){
for( int iz = zold + 1; iz <= znew; ++iz ){
  xv = xoffset;
  yv = yoffset;
  zv = zoffset;
  // z > 0 contribution
  xv += 1. * iz * t3[ 0 ];
  yv += 1. * iz * t3[ 1 ];
  zv += 1. * iz * t3[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  // z < 0 contribution
  xv -= 2. * iz * t3[ 0 ];
  yv -= 2. * iz * t3[ 1 ];
  zv -= 2. * iz * t3[ 2 ];
  xv2 = xv * xv;
  yv2 = yv * yv;
  zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
}
}
}
//
//
// now doing the (x,y,0) type translations, where x,y > 0
//
//
for( int ix = 1; ix <= xnew; ++ix ){
for( int iy = 1; iy <= ynew; ++iy ){
for( int iz = 0; iz <= 0; ++iz ){
  // making sure we're actually in the new box
  if( ix > xold || iy > yold || iz > zold ){
    x1 = 1. * ix * t1[ 0 ];
    x2 = 1. * iy * t2[ 0 ];
    x3 = 1. * iz * t3[ 0 ];

    y1 = 1. * ix * t1[ 1 ];
    y2 = 1. * iy * t2[ 1 ];
    y3 = 1. * iz * t3[ 1 ];

    z1 = 1. * ix * t1[ 2 ];
    z2 = 1. * iy * t2[ 2 ];
    z3 = 1. * iz * t3[ 2 ];

    xv = x3 + xoffset;
    yv = y3 + yoffset;
    zv = z3 + zoffset;
    // x > 0, y > 0 contribution
    xv += x1 + x2;
    yv += y1 + y2; 
    zv += z1 + z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y > 0 contribution
    xv -= 2. * x1; 
    yv -= 2. * y1; 
    zv -= 2. * z1;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y < 0 contribution
    xv -= 2. * x2; 
    yv -= 2. * y2; 
    zv -= 2. * z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y < 0 contribution
    xv += 2. * x1; 
    yv += 2. * y1; 
    zv += 2. * z1;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  }
}
}
}
for( int ix = 0; ix <= 0; ++ix ){
for( int iy = 1; iy <= ynew; ++iy ){
for( int iz = 1; iz <= znew; ++iz ){
  // making sure we're actually in the new box
  if( ix > xold || iy > yold || iz > zold ){
    x1 = 1. * ix * t1[ 0 ];
    x2 = 1. * iy * t2[ 0 ];
    x3 = 1. * iz * t3[ 0 ];

    y1 = 1. * ix * t1[ 1 ];
    y2 = 1. * iy * t2[ 1 ];
    y3 = 1. * iz * t3[ 1 ];

    z1 = 1. * ix * t1[ 2 ];
    z2 = 1. * iy * t2[ 2 ];
    z3 = 1. * iz * t3[ 2 ];

    xv = x1 + xoffset;
    yv = y1 + yoffset;
    zv = z1 + zoffset;
    // x > 0, y > 0 contribution
    xv += x3 + x2;
    yv += y3 + y2; 
    zv += z3 + z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y > 0 contribution
    xv -= 2. * x3; 
    yv -= 2. * y3; 
    zv -= 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y < 0 contribution
    xv -= 2. * x2; 
    yv -= 2. * y2; 
    zv -= 2. * z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y < 0 contribution
    xv += 2. * x3; 
    yv += 2. * y3; 
    zv += 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  }
}
}
}
for( int ix = 1; ix <= xnew; ++ix ){
for( int iy = 0; iy <= 0; ++iy ){
for( int iz = 1; iz <= znew; ++iz ){
  // making sure we're actually in the new box
  if( ix > xold || iy > yold || iz > zold ){
    x1 = 1. * ix * t1[ 0 ];
    x2 = 1. * iy * t2[ 0 ];
    x3 = 1. * iz * t3[ 0 ];

    y1 = 1. * ix * t1[ 1 ];
    y2 = 1. * iy * t2[ 1 ];
    y3 = 1. * iz * t3[ 1 ];

    z1 = 1. * ix * t1[ 2 ];
    z2 = 1. * iy * t2[ 2 ];
    z3 = 1. * iz * t3[ 2 ];

    xv = x2 + xoffset;
    yv = y2 + yoffset;
    zv = z2 + zoffset;
    // x > 0, y > 0 contribution
    xv += x3 + x1;
    yv += y3 + y1; 
    zv += z3 + z1;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y > 0 contribution
    xv -= 2. * x3; 
    yv -= 2. * y3; 
    zv -= 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y < 0 contribution
    xv -= 2. * x1; 
    yv -= 2. * y1; 
    zv -= 2. * z1;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y < 0 contribution
    xv += 2. * x3; 
    yv += 2. * y3; 
    zv += 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  }
}
}
}
//
//
// now doing the (x,y,z) type translations, where x,y,z > 0
//
//
for( int iz = 1; iz <= znew; ++iz ){
for( int ix = 1; ix <= xnew; ++ix ){
for( int iy = 1; iy <= ynew; ++iy ){
  // making sure we're actually in the new box
  if( ix > xold || iy > yold || iz > zold ){
    x1 = 1. * ix * t1[ 0 ];
    x2 = 1. * iy * t2[ 0 ];
    x3 = 1. * iz * t3[ 0 ];

    y1 = 1. * ix * t1[ 1 ];
    y2 = 1. * iy * t2[ 1 ];
    y3 = 1. * iz * t3[ 1 ];

    z1 = 1. * ix * t1[ 2 ];
    z2 = 1. * iy * t2[ 2 ];
    z3 = 1. * iz * t3[ 2 ];

    xv = xoffset;
    yv = yoffset;
    zv = zoffset;
    // x > 0, y > 0, z > 0 contribution
    xv += x1 + x2 + x3;
    yv += y1 + y2 + y3; 
    zv += z1 + z2 + z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y > 0, z < 0 contribution
    xv -= 2. * x3; 
    yv -= 2. * y3; 
    zv -= 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y < 0, z < 0 contribution
    xv -= 2. * x2; 
    yv -= 2. * y2; 
    zv -= 2. * z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x > 0, y < 0, z > 0 contribution
    xv += 2. * x3; 
    yv += 2. * y3; 
    zv += 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y < 0, z > 0 contribution
    xv -= 2. * x1; 
    yv -= 2. * y1; 
    zv -= 2. * z1;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y > 0, z > 0 contribution
    xv += 2. * x2; 
    yv += 2. * y2; 
    zv += 2. * z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y > 0, z < 0 contribution
    xv -= 2. * x3; 
    yv -= 2. * y3; 
    zv -= 2. * z3;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
    // x < 0, y < 0, z < 0 contribution
    xv -= 2. * x2; 
    yv -= 2. * y2; 
    zv -= 2. * z2;
    xv2 = xv * xv;
    yv2 = yv * yv;
    zv2 = zv * zv;
xyzvals[ counter * 3 + 0 ] = xv;
xyzvals[ counter * 3 + 1 ] = yv;
xyzvals[ counter * 3 + 2 ] = zv;
xyz_sq_vals[ counter * 3 + 0 ] = xv2;
xyz_sq_vals[ counter * 3 + 1 ] = yv2;
xyz_sq_vals[ counter * 3 + 2 ] = zv2;
counter++;
//printf( "<%6d> [ %3.0f, %3.0f, %3.0f ] \n", counter, xv, yv, zv );
  }
}
}
}

xold = xnew;
yold = ynew;
zold = znew;

//for( int i = 0; i < nvals; ++i ){
//  printf( "<%6d> [ %12.8f, %12.8f, %12.8f ] \n", i, xyzvals[ 3 * i + 0 ], xyzvals[ 3 * i + 1 ], xyzvals[ 3 * i + 2 ] );
//}

}

void generate_tvecs( 
  int& nvals,
  const Eigen::Matrix3d& Tmat,
  int in_xold, int in_xnew,
  std::vector< double >& xyzvals,
  std::vector< double >& xyz_sq_vals
){
  std::vector< double > tmatrix;
  tmatrix.resize( 9 );
  tmatrix[ 0 ] = Tmat( 0, 0 ); tmatrix[ 1 ] = Tmat( 0, 1 ); tmatrix[ 2 ] = Tmat( 0, 2 );
  tmatrix[ 3 ] = Tmat( 1, 0 ); tmatrix[ 4 ] = Tmat( 1, 1 ); tmatrix[ 5 ] = Tmat( 1, 2 );
  tmatrix[ 6 ] = Tmat( 2, 0 ); tmatrix[ 7 ] = Tmat( 2, 1 ); tmatrix[ 8 ] = Tmat( 2, 2 );
  generate_tvecs( nvals, tmatrix, in_xold, in_xnew, in_xold, in_xnew, in_xold, in_xnew,
                  0.0, 0.0, 0.0, xyzvals, xyz_sq_vals );
};
