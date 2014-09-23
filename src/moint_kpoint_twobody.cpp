#include <vector>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "cellinfo.h"
#include "common.h"
#include "moint_kpoint_twobody.h"
#include "Eigen/Dense"
#include <time.h>
using namespace std;
using namespace Eigen;

static vector< complex< double > > cos_sin_arg;

static int get_arg_index(
    int kpoint,
    int trans,
    int nkpt,
    int ntrans
){
    int index = trans + ntrans * kpoint;
    return index;
}

static void setup_cos_sin( 
    SuperCell& SCell
){
    int index;
    int nkpt = SCell.nkpt;
    int ntrans = SCell.total_number_of_cells;
    cos_sin_arg.resize( nkpt * ntrans );
    for( int ik = 0; ik < nkpt; ++ik ){
      for( int it = 0; it < ntrans; ++it ){
        index = get_arg_index( ik, it, nkpt, ntrans );
        real( cos_sin_arg[ index ] ) = cos( ( SCell.kpoints[ ik ] ).dot( SCell.translations[ it ] ) ); 
        imag( cos_sin_arg[ index ] ) = sin( ( SCell.kpoints[ ik ] ).dot( SCell.translations[ it ] ) ); 
      }
    }
}

std::vector< int > find_kpoint_conj( 
    std::vector< Eigen::Vector3i > reduced_kpoint,
    int nkx,
    int nky,
    int nkz,
    int nkpt
){
    vector< int > conj_loc;
    conj_loc.resize( nkpt );
    Vector3i kpt, ckpt;
    int k0, k1, k2;
    for( int i = 0; i < nkpt; ++i ){
      kpt = reduced_kpoint[ i ];
      for( int j = 0; j < nkpt; ++j ){
        ckpt = reduced_kpoint[ j ];
        k0 = abs( ckpt( 0 ) + kpt( 0 ) );
        if( k0 >= nkx ) k0 -= nkx;
        if( k0 != 0 ) break;
        k1 = abs( ckpt( 1 ) + kpt( 1 ) );
        if( k1 >= nky ) k1 -= nky;
        if( k1 != 0 ) break;
        k2 = abs( ckpt( 2 ) + kpt( 2 ) );
        if( k2 >= nkz ) k2 -= nkz;
        if( k2 != 0 ) break;
        conj_loc[ i ] = j; // since we have that ( kpt + ckpt == 0 )
      } 
    } 
    return conj_loc;
}


void moint_driver_k::kpoint_init(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXd& inevals,
    VMatrixXcd& inevecsXcd
){
    cout << "SETUP : KPOINT INTEGRAL DRIVER" << endl;
    evecsXcd = inevecsXcd;
    evals = inevals;
    nmo_ucell = UCell.nao;
    nkpt = SCell.nkpt;
    ntrans = SCell.total_number_of_cells;
    nmo_scell = SCell.nao;

    //atomic_dist_matr = inatomic_dist_matr;
    //ewald_pot_matr = inewald_pot_matr;
    kernel_matr.resize( nmo_scell, nmo_scell );
    for( int iao = 0; iao < nmo_ucell; ++iao ){
    for( int jao = 0; jao < nmo_ucell; ++jao ){
      for( int itrans = 0; itrans < ntrans; ++itrans ){
      for( int jtrans = 0; jtrans < ntrans; ++jtrans ){
        Vector3i tvec = SCell.reduced_t[ jtrans ] - SCell.reduced_t[ itrans ];
        bool found = false;
        int index = SCell.get_supercell_index( tvec( 0 ), tvec( 1 ), tvec( 2 ), &found );
        kernel_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) = aoints.getXCInt( iao, jao, index );
      }
      }
    }
    }

    cout << " - initializing cos and sin arguments..." << flush;
    setup_cos_sin( SCell ); 
    cout << "done." << endl;
    cout << "SETUP COMPLETE!" << endl;
}

void moint_driver_k::write_out_moints( std::vector< std::complex< double > > kmoints, const char* ofilename ){
    ofstream outfile( ofilename, std::ios::binary );
    size_t kmo_size = kmoints.size();
    outfile.write( (char*) &kmo_size, sizeof( size_t ) );
    outfile.write( (char*) &kmoints, kmoints.size() * sizeof( complex< double > ) );
    outfile.close();
}

std::vector< std::complex< double > > moint_driver_k::read_in_moints( const char* ifilename ){
    std::vector< std::complex< double > > kmoints;
    ifstream infile( ifilename, std::ios::binary );
    size_t kmo_size;
    infile.read( (char*) &kmo_size, sizeof( size_t ) );
    kmoints.resize( kmo_size );
    infile.read( (char*) &kmoints, kmoints.size() * sizeof( complex< double > ) );
    infile.close();
    return kmoints;
}

// Only have integrals of the form (G_i G_j | G_k G_l) where [ G_i - G_k + G_j - G_l ] =equiv= 0 ].
// From i,j,k specifications there can thus only be one G_l, so only need to specify i,j,k for kpoints
// Also have the symmetry < I J | K L > = < J* I* | L* K* >, but let's ignore that for now...

int get_ijkl_kpoint(
    int kindex,
    int orbindex,
    int max_k,
    int max_orb
){
    int index = kindex + max_k * orbindex;
    return index;
}

void moint_driver_k::check_exchange_energy( 
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    vector< complex< double > >& moints 
){
    int moindex, nindex_k, nindex_ijkl;
    nindex_k = nkpt * nkpt * nkpt;
    nindex_ijkl = nmo_ucell * nmo_ucell * nmo_ucell * nmo_ucell; 
    printf( "CHECKING EXCHANGE ENERGY USING MO INTEGRALS...\n" );
    complex< double > xc;
    xc *= 0.0;
    int ind1, ind2;
    for( int i = 0; i < nkpt; ++i ){
      for( int j = 0; j < nkpt; ++j ){
        for( int iao = 0; iao < (int)( nmo_ucell / 2 ); ++iao ){
          for( int jao = 0; jao < (int)( nmo_ucell / 2 ); ++jao ){
            ind1 = i + j * nkpt + j * nkpt * nkpt;
            ind2 = iao + jao * nmo_ucell + jao * nmo_ucell * nmo_ucell + iao * nmo_ucell * nmo_ucell * nmo_ucell;
            moindex = get_ijkl_kpoint( ind1, ind2, nindex_k, nindex_ijkl );
            cout << "GETTING INDEX : " << moindex << " " << moints[ moindex ] << endl;
            xc += - 1. * moints[ moindex ];
          }
        }
      }
    }
    printf( "  - EXCHANGE ENERGY        (re,im): %20.16f %20.16f \n", real( xc ), imag( xc ) );
    xc /= nkpt;
    printf( "  - EXCHANGE ENERGY / nkpt (re,im): %20.16f %20.16f \n", real( xc ), imag( xc ) );
}

vector< complex< double > > moint_driver_k::create_kpoint_2e_ints(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXcd& inevecsXcd
){
    int index_k, nindex_k; 
    int index_ijkl, nindex_ijkl; 
    int moindex;
    bool use_contr;
    vector< complex< double > > outvec;
    nindex_k = nkpt * nkpt * nkpt;
    nindex_ijkl = nmo_ucell * nmo_ucell * nmo_ucell * nmo_ucell; 
    int mosize = nindex_ijkl * nindex_k;
    cout << "creating mo integrals with size : " << mosize << endl;
    outvec.resize( mosize );
    for( int i = 0; i < mosize; ++i ) outvec[ i ] *= 0.0;
    complex< double > c1, c2, c3, c4;
    complex< double > exp1, exp2, exp3, exp4;
    complex< double > sum1, sum2, sum3, sum4;
    clock_t starttime = clock();

    bool found;
    int iarg;
    int lkp;
    int nmo_ucell_sq = nmo_ucell * nmo_ucell;
    int nmo_ucell_cu = nmo_ucell * nmo_ucell * nmo_ucell;
    // ... these are the k-indices
    for( int ikp = 0; ikp < nkpt; ++ikp ){ 
    for( int jkp = 0; jkp < nkpt; ++jkp ){ 
    for( int kkp = 0; kkp < nkpt; ++kkp ){ 
      index_k = ikp + jkp * nkpt + kkp * nkpt * nkpt;
      // ... G_l = G_i - G_j + G_k
      Vector3i lkpoint = SCell.reduced_k[ ikp ] - SCell.reduced_k[ jkp ] + SCell.reduced_k[ kkp ];

      if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
      if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
      if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;

      if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
      if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
      if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
      
      lkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );

      //cout << "K_1 : " << SCell.reduced_k[ ikp ]( 0 ) << ", " << SCell.reduced_k[ ikp ]( 1 ) << ", " << SCell.reduced_k[ ikp ]( 2 ) << endl;
      //cout << "K_3 : " << SCell.reduced_k[ kkp ]( 0 ) << ", " << SCell.reduced_k[ kkp ]( 1 ) << ", " << SCell.reduced_k[ kkp ]( 2 ) << endl;
      //cout << "K_2 : " << SCell.reduced_k[ jkp ]( 0 ) << ", " << SCell.reduced_k[ jkp ]( 1 ) << ", " << SCell.reduced_k[ jkp ]( 2 ) << endl;
      //cout << "K_4 : " << SCell.reduced_k[ lkp ]( 0 ) << ", " << SCell.reduced_k[ lkp ]( 1 ) << ", " << SCell.reduced_k[ lkp ]( 2 ) << endl;

      //cout << "FOR K POINTS : " << ikp << ", " << jkp << ", " << kkp << ", " << lkp << " : " << endl;
      if( !found ){
        cout << "KPOINT : " << lkpoint( 0 ) << ", " << lkpoint( 1 ) << ", " << lkpoint( 2 ) << " not found???" << endl;
        exit( EXIT_FAILURE );
      }
      // ... at this k-point, we have nmo_ucell MO's that are linear combinations of bloch-wavefunctions
      for( int imo = 0; imo < nmo_ucell; ++imo ){
      for( int jmo = 0; jmo < nmo_ucell; ++jmo ){
      for( int kmo = 0; kmo < nmo_ucell; ++kmo ){
      for( int lmo = 0; lmo < nmo_ucell; ++lmo ){
        index_ijkl = imo + jmo * nmo_ucell + kmo * nmo_ucell_sq + lmo * nmo_ucell_cu;
        //cout << "CREATING MO's FOR : " << imo << ", " << jmo << ", " << kmo << ", " << lmo << " : " << index_ijkl << endl;
        moindex = get_ijkl_kpoint( index_k, index_ijkl, nindex_k, nindex_ijkl );
        //cout << "index : " << moindex << endl;
        // ... now we have our sum over the ao's in the unit cell
        // ... use PPP properties so that only 2 summations necessary
        sum1 *= 0.0;
        for( int iao = 0; iao < nmo_ucell; ++iao ){  
          c1 = evecsXcd.irrep( ikp )( iao, imo );
          c2 = evecsXcd.irrep( jkp )( iao, jmo );
          sum2 *= 0.0;
          for( int jao = 0; jao < nmo_ucell; ++jao ){  
            c3 = evecsXcd.irrep( kkp )( jao, kmo );
            c4 = evecsXcd.irrep( lkp )( jao, lmo );
            // ... now summing over translational vectors, again using the PPP properties 
            sum3 *= 0.0;
            //cout << "coefficients : " << c1 << ", " << c2 << ", " << c3 << ", " << c4 << endl;
            for( int itv = 0; itv < ntrans; ++itv ){
              iarg = get_arg_index( ikp, itv, nkpt, ntrans );
              exp1 = cos_sin_arg[ iarg ]; 
              iarg = get_arg_index( jkp, itv, nkpt, ntrans );
              exp2 = cos_sin_arg[ iarg ]; 
              sum4 *= 0.0;
              for( int jtv = 0; jtv < ntrans; ++jtv ){
                iarg = get_arg_index( kkp, jtv, nkpt, ntrans ); 
                exp3 = cos_sin_arg[ iarg ]; 
                iarg = get_arg_index( lkp, jtv, nkpt, ntrans ); 
                exp4 = cos_sin_arg[ iarg ]; 
                sum4 += conj( exp3 ) * exp4 * kernel_matr( iao + itv * nmo_ucell, jao + jtv * nmo_ucell );
                // ... we are onsite!
                if( iao == jao && itv == jtv ) 
                  sum4 += conj( exp3 ) * exp4 * aoints.getHubbardU( UCell.ElementStr[ iao ] );
//                cout << " sum4 : " << sum4 << endl;
              }
              sum3 += conj( exp1 ) * exp2 * sum4;
            }
            sum2 += conj( c3 ) * c4 * sum3;
          }
          sum1 += conj( c1 ) * c2 * sum2;
        }
        sum1 /= ( 1. * nkpt * nkpt );
        //cout << "VAL : " << sum1 << endl;
        outvec[ moindex ] = sum1;
      }
      }
      }
      }
    }
    }
    }

   cout << "elapsed time : " << double( clock() - starttime ) / (double)CLOCKS_PER_SEC << " seconds " << endl;

    return outvec;
}

//
//
// For a good idea of what's going on.. check the 'create_kpoint_2e_ints' driver.. Much more organized
// than this one.  Basically the only change is that this fast one does intermediate contractions.
//
//
vector< complex< double > > moint_driver_k::create_kpoint_2e_ints_fast_FULL(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXcd& inevecsXcd
){
    vector< complex< double > > integrals;
    integrals = create_kpoint_2e_ints_fast( UCell, SCell, aoints, inevecsXcd, 0, UCell.nao, 0, UCell.nao, 0, UCell.nao, 0, UCell.nao );
    return integrals;
}

vector< complex< double > > moint_driver_k::create_kpoint_2e_ints_fast(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    VMatrixXcd& inevecsXcd,
    int pstart,
    int plast,
    int qstart,
    int qlast,
    int rstart,
    int rlast,
    int sstart,
    int slast
){
    int index_k, nindex_k; 
    int index_ijkl, nindex_ijkl; 
    int moindex;
    bool use_contr;
    vector< complex< double > > outvec;
    nindex_k = nkpt * nkpt * nkpt;
    nindex_ijkl = nmo_ucell * nmo_ucell * nmo_ucell * nmo_ucell; 
    int mosize = nindex_ijkl * nindex_k;
    cout << "creating mo integrals with size : " << mosize << endl;
    outvec.resize( mosize );
    for( int i = 0; i < mosize; ++i ) outvec[ i ] *= 0.0;
    complex< double > c1, c2, c3, c4;
    complex< double > exp1, exp2, exp3, exp4;
    complex< double > sum1, sum2, sum3, sum4;
    clock_t starttime = clock();

    int counter = 0;
    int counter2 = 0;
    vector< complex< double > > store_ij;
    complex< double > tempstore1, tempstore2;
    store_ij.resize( nmo_ucell * ntrans );
    bool found;
    int iarg;
    int lkp;
    int nmo_ucell_sq = nmo_ucell * nmo_ucell;
    int nmo_ucell_cu = nmo_ucell * nmo_ucell * nmo_ucell;
    // ... these are the k-indices
    for( int ikp = 0; ikp < nkpt; ++ikp ){ 
    for( int jkp = 0; jkp < nkpt; ++jkp ){ 
      printf( "CREATING ALL MO's FOR K_i K_j : %4d %4d \n", ikp, jkp );
      // ... at this k-point, we have nmo_ucell MO's that are linear combinations of bloch-wavefunctions
      for( int imo = pstart; imo < plast; ++imo ){
      for( int jmo = qstart; jmo < qlast; ++jmo ){
      // let's construct all the values for this IJ
        counter = 0;
        for( int i = 0; i < store_ij.size(); ++i ) 
          store_ij[ i ] *= 0.0;
        for( int itv = 0; itv < ntrans; ++itv ){
          iarg = get_arg_index( ikp, itv, nkpt, ntrans );
          exp1 = cos_sin_arg[ iarg ]; 
          iarg = get_arg_index( jkp, itv, nkpt, ntrans );
          exp2 = cos_sin_arg[ iarg ]; 
          tempstore1 = conj( exp1 ) * exp2;
          for( int iao = 0; iao < nmo_ucell; ++iao, ++counter ){  
            c1 = evecsXcd.irrep( ikp )( iao, imo );
            c2 = evecsXcd.irrep( jkp )( iao, jmo );
            tempstore2 = conj( c1 ) * c2 * tempstore1;
            counter2 = 0;
            for( int jtv = 0; jtv < ntrans; ++jtv ){
              for( int jao = 0; jao < nmo_ucell; ++jao, ++counter2 ){
                store_ij[ counter2 ] += tempstore2 * kernel_matr( counter, counter2 );
                if( iao == jao && itv == jtv ) // we are onsite
                    store_ij[ counter2 ] += tempstore2 * aoints.getHubbardU( UCell.ElementStr[ iao ] );
              }
            }
          }
        }


      for( int kkp = 0; kkp < nkpt; ++kkp ){ 
        index_k = ikp + jkp * nkpt + kkp * nkpt * nkpt;
        // ... G_l = G_i - G_j + G_k
        Vector3i lkpoint = SCell.reduced_k[ ikp ] - SCell.reduced_k[ jkp ] + SCell.reduced_k[ kkp ];
     
        if( lkpoint( 0 ) <  0         ) lkpoint( 0 ) += SCell.nkx;
        if( lkpoint( 1 ) <  0         ) lkpoint( 1 ) += SCell.nky;
        if( lkpoint( 2 ) <  0         ) lkpoint( 2 ) += SCell.nkz;
     
        if( lkpoint( 0 ) >= SCell.nkx ) lkpoint( 0 ) -= SCell.nkx;
        if( lkpoint( 1 ) >= SCell.nky ) lkpoint( 1 ) -= SCell.nky;
        if( lkpoint( 2 ) >= SCell.nkz ) lkpoint( 2 ) -= SCell.nkz;
        
        lkp = SCell.get_kpoint_index( lkpoint( 0 ), lkpoint( 1 ), lkpoint( 2 ), &found );
     
        if( !found ){
          cout << "KPOINT : " << lkpoint( 0 ) << ", " << lkpoint( 1 ) << ", " << lkpoint( 2 ) << " not found???" << endl;
          exit( EXIT_FAILURE );
        }


        for( int kmo = rstart; kmo < rlast; ++kmo ){
        for( int lmo = sstart; lmo < slast; ++lmo ){
          index_ijkl = imo + jmo * nmo_ucell + kmo * nmo_ucell_sq + lmo * nmo_ucell_cu;
          moindex = get_ijkl_kpoint( index_k, index_ijkl, nindex_k, nindex_ijkl );
          counter = 0;
          sum3 *= 0.0;
          for( int jtv = 0; jtv < ntrans; ++jtv ){
              iarg = get_arg_index( kkp, jtv, nkpt, ntrans ); 
              exp3 = cos_sin_arg[ iarg ]; 
              iarg = get_arg_index( lkp, jtv, nkpt, ntrans ); 
              exp4 = cos_sin_arg[ iarg ]; 
              sum4 *= 0.0;
              for( int jao = 0; jao < nmo_ucell; ++jao, ++counter ){  
                c3 = evecsXcd.irrep( kkp )( jao, kmo );
                c4 = evecsXcd.irrep( lkp )( jao, lmo );
                sum4 += conj( c3 ) * c4 * store_ij[ counter ];
              }
              sum3 += conj( exp3 ) * exp4 * sum4;
          }
          sum3 /= ( 1. * nkpt * nkpt );
          //cout << "VAL : " << sum1 << endl;
          outvec[ moindex ] = sum3;
        }
        }



      }
      }
    }
    }
    }

   cout << "elapsed time : " << double( clock() - starttime ) / (double)CLOCKS_PER_SEC << " seconds " << endl;

    return outvec;
}
