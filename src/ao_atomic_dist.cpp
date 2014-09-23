#include "ao_ints.h"
#include "cellinfo.h"
#include "stdio.h"
using namespace std;
using namespace Eigen;


void aoIntegralFactory::makeDistMatr(
    UnitCell& UCell,
    SuperCell& SCell,
    bool explicit_MIC
){
    V3d distvec;

    if( !explicit_MIC )
        printf( "Warning : makeDistMatr does not have explicit_MIC on, may not find minimum image! \n" );
 
    Vector3d v1;
    RowVector3d vtemp;
    vector< double > frac_coords;
    if( dim > 0 )
        frac_coords.resize( dim );
    double vproj, dist, distsq, tempdsq, tempd;
    Vector3d vtranslate;
    MatrixXd vtranspose;
    vtranspose.resize( 1, 3 );
    int distances_changed = 0;
    bool smallerfound;


    /* Although we may have translation vectors defined, the wigner-seitz cell may not be formed
       by a translation vector in each direction (for very skewed lattices).  The following sets
       up some translation vectors to search over in order to find the minimum image */

    const int t1max = 3;
    const int t2max = 3;
    const int t3max = 3;
    vector< RowVector3d > ptrans;
    ptrans.resize( ( 2 * t1max + 1 ) * ( 2 * t2max + 1 ) * ( 2 * t3max + 1 ) );
    int vcount = 0;
    vector< double > ptsq;
    vector< double > pt;
    ptsq.resize( 3 * ptrans.size() );
    pt.resize( 3 * ptrans.size() );
    for( int i = -t1max; i <= t1max; ++i ){
      for( int j = -t2max; j <= t2max; ++j ){
        for( int k = -t3max; k <= t3max; ++k, ++vcount ){
          vtemp = 1. * i * SCell.T.row( 0 ) + 1. * j * SCell.T.row( 1 ) + 1. * k * SCell.T.row( 2 );
          ptrans[ vcount ] = vtemp;
        }
      }
    }

    /* now we search for the nearest images */

    int naoUsqp1 = ( naoUnitCell * ( naoUnitCell + 1 ) ) / 2;
    if( aoDistMatr.empty() )
      aoDistMatr.resize( naoUsqp1 * nTrans );

    for( int i = 0; i < naoUnitCell; i++ ){
       for( int itrans = 0; itrans < nTrans; itrans++ ){
       for( int j = i; j < naoUnitCell; j++ ){
           smallerfound = false;
           vtranslate = SCell.translations[ itrans ] + UCell.coords[ j ] - UCell.coords[ i ]; 
           vtranspose = vtranslate.transpose();
           v1 = vtranslate;
           for( int ii = 0; ii < dim; ii++ ){
              frac_coords[ ii ] = invmat( ii, 0 ) * v1( 0 ) + invmat( ii, 1 ) * v1( 1 ) + invmat( ii, 2 ) * v1( 2 );
              if( frac_coords[ ii ] < 0.0 )
                v1 -= round( frac_coords[ ii ] + 1e-15 * frac_coords[ ii ] ) * SCell.T.row( ii );
              else
                v1 -= round( frac_coords[ ii ] - 1e-15 * frac_coords[ ii ] ) * SCell.T.row( ii );
           }
           distsq = v1.dot( v1 );
           dist = sqrt( distsq );

           if( explicit_MIC ){
                for( int i = 0; i < ptrans.size(); ++i ){
                    vtemp = vtranspose + ptrans[ i ]; 
                    tempdsq = vtemp.dot( vtemp );
                    tempd = sqrt( tempdsq );
                    if( tempd < ( dist - 1e-5 ) ){
                        dist = tempd;
                        smallerfound = true; 
                    }
                }
                if( smallerfound ) distances_changed++;
           }

         aoDistMatr[ getPerElement( i, j, itrans ) ] = dist;
         //printf( "%20.16f  iat %3d  jat %3d  trans %3d    INDEX : %3d \n", dist, i, j, itrans, getPerElement( i, j, itrans ) );
//           aoDistMatr( i, j ) = dist;
//           aoDistMatr( j, i ) = dist;
       }
       }
    }

//    for( int i = 0; i < SCell.total_number_of_cells; ++i ){
//        for( int j = i; j < SCell.total_number_of_cells; ++j ){
//            aoDistMatr.block( i * nmo_ucell, j * nmo_ucell, nmo_ucell, nmo_ucell ) = \
//            aoDistMatr.block( 0, scell_matr( i, j ) * nmo_ucell, nmo_ucell, nmo_ucell );
//        }
//        for( int j = 0; j < i; ++j ){
//            aoDistMatr.block( i * nmo_ucell, j * nmo_ucell, nmo_ucell, nmo_ucell ) = \
//            aoDistMatr((*atomic_dist_matr).block( j * nmo_ucell, i * nmo_ucell, nmo_ucell, nmo_ucell ) ).transpose();
//        }
//    }

//    cout << "NUMBER OF DISTANCES DIFFERENT FROM 'MIC' METHOD : " << distances_changed << " ( out of " << ( nmo_scell * ( nmo_scell + 1 )) / 2 << " ) " << endl;
//    cout << "ATOMIC DISTANCE MATRIX : " << endl;
//    cout << (*atomic_dist_matr) << endl << endl;
}
