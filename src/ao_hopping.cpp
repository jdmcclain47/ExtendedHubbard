#include "ao_ints.h"
#include "cellinfo.h"
#include "ppp.h"
#include "math.h"
#include <cassert>
#include <cstdio> 
#include <vector>

void aoIntegralFactory::makeHoppingMatr(
  UnitCell& UCell,
  SuperCell& SCell
){
  assert( aoDistMatr.size() > 0 && "When calling setPPPHopptinMatr, aoDistMatr should already be set up!!!" );
  if( aoHopMatr.empty() )
    aoHopMatr.resize( aoDistMatr.size() ); 
  double dist;
  int index;
  for( int i = 0; i < naoUnitCell; ++i ){
  for( int j = i; j < naoUnitCell; ++j ){
    for( int itrans = 0; itrans < nTrans; ++itrans ){
      index = getPerElement( i, j, itrans );
      aoHopMatr[ index ] = PPP.get_hopping( UCell.ElementStr[ i ], UCell.ElementStr[ j ], aoDistMatr[ index ] );
      //printf( "ELEMENT_I: %3s  ELEMENT_J: %3s  DISTANCE: %12.8f  VALUE: %12.8f \n", UCell.ElementStr[ i ].c_str(), 
      //        UCell.ElementStr[ j ].c_str(), aoDistMatr[ index ], aoHopMatr[ index ] );
    } 
  }
  }
}
