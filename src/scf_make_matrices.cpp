#include <iostream>
#include "cellinfo.h"
#include "common.h"
#include "ao_ints.h"
#include "scf.h"
#include "Eigen/Dense"
using namespace std;

Eigen::MatrixXd SCF::make_coulomb_matr( 
  UnitCell& UCell, SuperCell& SCell, 
  aoIntegralFactory& aoints, VMatrixXd& inDens, 
  int unit_cell_image 
){
  Eigen::MatrixXd outmatr;
  outmatr = Eigen::MatrixXd::Zero( nmo_ucell, nmo_ucell );
  if( unit_cell_image == 0 ){
    for( int i = 0; i < nmo_ucell; ++i ){
      //
      //
      // note that when we are on-site, we need 4 different contributions :
      //         1)  the interaction with the on-site electron
      //         2)  the interaction with the on-site electron's replicas
      //         3)  the interaction with the on-site nucleus 
      //         4)  the interaction with the on-site nucleus's replicas
      // The (1) is taken care of here and (2) and (4) are taken care of in the next step.  This leave us
      // with (3) to be taken care of, but since it's just a constant it can absorb into the OnsiteE terms.
      // Note that it is done this way just to make a clear connection to the hubbard model and to be able
      // to switch back between the hubbard model with ease.
      //
      //
      outmatr( i, i ) += inDens.irrep( 0 )( i, i ) * aoints.getHubbardU( UCell.ElementStr[ i ] );
      for( int j = 0; j < nmo_ucell; ++j ){
        // 
        //
        // note that we still "add" the contribution when i==j and Rn==0, i.e. we are on the same site,
        // however, this is adding the potential due to its own background of replicated charges
        // for ewald type methods.  For cut-off type schemes, the getCoulombInt returns 0, and so 
        // the hubbard U is only counted once
        // 
        //
        for( int icell = 0; icell < SCell.total_number_of_cells; ++icell ){
          outmatr( i, i ) += ( inDens.irrep( 0 )( j, j ) - UCell.charge[ j ] ) * aoints.getCoulombInt( i, j, icell );
        }
      }
    }
  }else{
  // this is 0 according to the PPP model
  }
  return outmatr;
}


Eigen::MatrixXd SCF::make_exchange_matr( 
  UnitCell& UCell, SuperCell& SCell, 
  aoIntegralFactory& aoints, VMatrixXd& inDens, 
  int unit_cell_image 
){
  double nonehalf = -0.5;
  Eigen::MatrixXd outmatr;
  outmatr = Eigen::MatrixXd::Zero( nmo_ucell, nmo_ucell );
  for( int i = 0; i < nmo_ucell; ++i ){
    if( unit_cell_image == 0 )
      outmatr( i, i ) += inDens.irrep( 0 )( i, i ) * aoints.getHubbardU( UCell.ElementStr[ i ] );
    for( int j = 0; j < nmo_ucell; ++j ){
      // note that we still "add" the contribution when i==j and Rn==0, i.e. we are on the same site,
      // however, this is adding the potential due to its own background of replicated charges
      // for ewald type methods.  For cut-off type schemes, the getXCInt returns 0, and so 
      // the hubbard U is only counted once. For a further explanation see make_coulomb_matr.
      outmatr( i, j ) += inDens.irrep( unit_cell_image )( i, j ) * aoints.getXCInt( i, j, unit_cell_image );
    }
  }
  outmatr *= nonehalf;
  return outmatr;
}


Eigen::MatrixXd SCF::make_kinetic_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, int unit_cell_image ){
  Eigen::MatrixXd outmatr;
  outmatr = Eigen::MatrixXd::Zero( nmo_ucell, nmo_ucell );
  for( int i = 0; i < nmo_ucell; ++i ){
    for( int j = 0; j < nmo_ucell; ++j ){
      outmatr( i, j ) += aoints.getHopping( i, j, unit_cell_image );
    }
  }
  return outmatr;
}


Eigen::MatrixXd SCF::make_onsiteE_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, int unit_cell_image ){
  Eigen::MatrixXd outmatr;
  outmatr = Eigen::MatrixXd::Zero( nmo_ucell, nmo_ucell );
  if( unit_cell_image == 0 ){
    for( int i = 0; i < nmo_ucell; ++i ){
      outmatr( i, i ) += aoints.getOnsiteE( UCell.ElementStr[ i ] );
    }
  }else{
    // do nothing, since onsiteE only serves to shift the diagonal
  }
  return outmatr;
}
