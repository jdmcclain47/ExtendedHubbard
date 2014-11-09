#include "ao_ints.h"
#include "cellinfo.h"
#include "ppp.h"
#include "math.h"
#include <cassert>
#include <cstdio> 
#include <vector>

void aoIntegralFactory::makeCutMatrWS(
  UnitCell& UCell,
  SuperCell& SCell
){
  assert( aoDistMatr.size() > 0 && "When calling makeCutMatr, aoDistMatr should already be set up!!!" );
  if( aoCutMatrWS.empty() )
    aoCutMatrWS.resize( aoDistMatr.size() ); 
  double dist;
  int index;
  for( int i = 0; i < naoUnitCell; ++i ){
  for( int j = i; j < naoUnitCell; ++j ){
    for( int itrans = 0; itrans < nTrans; ++itrans ){
      index = getPerElement( i, j, itrans );
      if( i == j && itrans == 0 ) // then we are on the same site, and so don't add anything
        aoCutMatrWS[ index ] = 0.0; 
      else
        aoCutMatrWS[ index ] = 1./aoDistMatr[index]; //getPPPInt( UCell.ElementStr[ i ], UCell.ElementStr[ j ], aoDistMatr[ index ] );
      //printf( "ELEMENT_I: %3s  ELEMENT_J: %3s  DISTANCE: %12.8f  VALUE: %12.8f \n", UCell.ElementStr[ i ].c_str(), 
      //        UCell.ElementStr[ j ].c_str(), aoDistMatr[ index ], aoCutMatrWS[ index ] );
    } 
  }
  }
}


