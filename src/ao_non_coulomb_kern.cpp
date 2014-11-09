#include "Eigen/Dense"
#include "ao_ints.h"
#include "cellinfo.h"
#include "find_bounding_sphere.h"
#include "ppp.h"
#include "math.h"
#include <cassert>
#include <cstdio> 
#include <vector>

void aoIntegralFactory::makeNonCoulombMatr(
  UnitCell& UCell,
  SuperCell& SCell,
  double rcut
){
  assert( aoDistMatr.size() > 0 && "When calling setNonCoulombMatr, aoDistMatr should already be set up!!!" );
  if( pppkern == COULOMB )
    return; // no need to create a correction to the coulomb kernel
  if( aoNonCoulombMatr.empty() )
    aoNonCoulombMatr.resize( aoDistMatr.size() ); 
  int nx, ny, nz;
  Eigen::Matrix3d tempTmat = SCell.T;
  findBoundingSphere( tempTmat, rcut, nx, ny, nz );
  nx++;
  ny++;
  nz++;
  printf( "MAKING CORRECTION TO COULOMB OPERATOR...\n" );
  printf( "  - correcting operator up to a distance : %20.16f \n", rcut );
  printf( "  - search space for trans. vectors      : (%3d,%3d,%3d) \n", nx, ny, nz );

  Eigen::RowVector3d tvec;
  Eigen::Vector3d dvec;
  double dist;
  int index;
  for( int i = 0; i < naoUnitCell; ++i ){
  for( int j = i; j < naoUnitCell; ++j ){
    for( int itrans = 0; itrans < nTrans; ++itrans ){
      index = getPerElement( i, j, itrans );
      aoNonCoulombMatr[ index ] = 0.0; 
      dvec = ( UCell.coords[ j ] + SCell.translations[ itrans ] - UCell.coords[ i ] );
      for( int inx = -nx; inx <= nx; ++inx ){
      for( int iny = -ny; iny <= ny; ++iny ){
      for( int inz = -nz; inz <= nz; ++inz ){
        tvec = 1. * inx * tempTmat.row( 0 ) + 1. * iny * tempTmat.row( 1 ) + 1. * inz * tempTmat.row( 2 );
        tvec( 0 ) += dvec( 0 );
        tvec( 1 ) += dvec( 1 );
        tvec( 2 ) += dvec( 2 );
        dist = sqrt( tvec.dot( tvec ) );
        /*
          now we make sure we are less than the rcut value and that we are not onsite
          since onsite/U values are taken care elsewhere ...
        */
        if( dist <= rcut && ( abs( itrans ) + abs( inx ) + abs( iny ) + abs( inz ) + abs( j - i ) != 0 ) ){
          aoNonCoulombMatr[ index ] += getPPPInt( UCell.ElementStr[ i ], UCell.ElementStr[ j ], dist ) - \
                                       1. / dist; 
        }
      } 
      } 
      } 
    } 
  }
  }
}
