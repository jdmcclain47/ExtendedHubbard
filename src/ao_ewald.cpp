#include "ewald.h"
#include <vector>
#include <math.h>
#include "cellinfo.h"
#include "ao_ints.h"
#include "Eigen/Dense"
using namespace Eigen;

void aoIntegralFactory::makeEwaldMatr(
    UnitCell& UCell,
    SuperCell& SCell,
    int tolEwald
){

  if( aoEwaldMatr.empty() )
    aoEwaldMatr.resize( naoUnit_sqp1 * nTrans );
  
  int index;
  Eigen::Vector3d coordj;
  std::vector< double > charge ( 1, 1.0 );
  for( int i = 0; i < naoUnitCell; ++i ){
  for( int j = i; j < naoUnitCell; ++j ){
    for( int itrans = 0; itrans < nTrans; ++itrans ){
      coordj = SCell.translations[ itrans ] + UCell.coords[ j ];
      index = getPerElement( i, j, itrans );
      //aoEwaldMatr[ index ] = potential_ewald_converge( UCell.coords[ i ], SCell.T, SCell.K, coordj, charge, tolEwald );
      if( i != j ){
        aoEwaldMatr[ index ] = potential_ewald_converge_FAST( UCell.coords[ i ], SCell.T, SCell.K, coordj, charge, tolEwald );
      }else{
        if( i > 0 ){
          aoEwaldMatr[ index ] = aoEwaldMatr[ getPerElement( 0, 0, itrans ) ]; 
        }else{
          aoEwaldMatr[ index ] = potential_ewald_converge_FAST( UCell.coords[ i ], SCell.T, SCell.K, coordj, charge, tolEwald );
        }
      }
    } 
  }
  }

}
