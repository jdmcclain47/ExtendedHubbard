#include <fstream>
#include <vector>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "cellinfo.h"
#include "common.h"
#include "mointegrals.h"
#include "Eigen/Dense"
#include <time.h>
using namespace std;
using namespace Eigen;

static size_t get_ijkl_index( int imo, int jmo, int kmo, int lmo ){
    size_t ij = ( imo > jmo ) ? ( imo * ( imo + 1 ) ) / 2 + jmo
                              : ( jmo * ( jmo + 1 ) ) / 2 + imo;
    size_t kl = ( kmo > lmo ) ? ( kmo * ( kmo + 1 ) ) / 2 + lmo
                              : ( lmo * ( lmo + 1 ) ) / 2 + kmo;
    size_t index = ( ij > kl ) ? ( ij * ( ij + 1 ) ) / 2 + kl
                            : ( kl * ( kl + 1 ) ) / 2 + ij;
    return index;
}

std::vector< double > moIntegralFactory::read_gamma_moints_from_file( const char* ifilename ){
    int in_nmo_scell;
    double val;
    int index;
    ifstream infile( ifilename, std::ios::binary );
    infile.read( (char*) &in_nmo_scell, sizeof( int ) );

    vector< double > outvec;
    size_t outsize = ( in_nmo_scell * ( in_nmo_scell - 1 ) ) / 2 + ( in_nmo_scell - 1 );
    outsize = ( outsize * ( outsize + 1 ) ) / 2 + outsize + 1;
    outvec.resize( outsize );
    for( int i = 0; i < outsize; ++i ) outvec[ i ] = 0.0;

    while( !infile.eof() ){
      infile.read( (char*) &index, sizeof( index ) );
      infile.read( (char*) &val, sizeof( double ) );
      outvec[ index ] = val;
    } 
    return outvec; 
}

double moIntegralFactory::make_2eint_pqrs( 
    const Eigen::MatrixXd& kernel_matr,
    const Eigen::MatrixXd& evecsXd,
    int p, int q, int r, int s
){
    double val;
    Eigen::VectorXd v1, v2, v3;
    v1.resize( nmo_scell );
    v2.resize( nmo_scell );
    v1 = evecsXd.col( p );
    v2 = evecsXd.col( q );
    for( int i = 0; i < nmo_scell; ++i ) v1( i ) = v1( i ) * v2( i );
    v2 = evecsXd.col( r );
    v3 = evecsXd.col( s );
    for( int i = 0; i < nmo_scell; ++i ) v2( i ) = v2( i ) * v3( i );
    val = v1.transpose() * kernel_matr * v2;
    return val;
}

void moIntegralFactory::write_gamma_moint_to_file( const char* ofilename ){
    int ij_index, kl_index, mo_index;
    double val;
    ofstream outfile( ofilename, std::ios::binary );
    //double dtol = pow( 10., 1. * itol );
    outfile.write( (char*) &nmo_scell, sizeof( int ) );


    bool create_new_kl;
    bool create_first_cont, create_second_cont, create_third_cont, create_fourth_cont;
    for( int i = 0; i < nmo_scell; ++i ){
    create_first_cont = true;
    printf( "Performing gamma transform ; Completed (%3d /%3d ) ... \n", (i+1), nmo_scell ); cout << flush;
    for( int j = 0; j <= i; ++j ){
      create_second_cont = true;
      ij_index = ( i * ( i + 1 ) ) / 2 + j;
      create_new_kl = true;
      for( int k = 0; k < nmo_scell; ++k ){
      create_third_cont = true;
      for( int l = 0; l <= k; ++l ){
        create_fourth_cont = true;
        kl_index = ( k * ( k + 1 ) ) / 2 + l;
//        if( ij_index >= kl_index ){
        if( i >= k ){
          mo_index = ( ij_index * ( ij_index + 1 ) ) / 2 + kl_index;
          val = quarter_gamma_transform( i, j, k, l, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
          //outfile.write( (char*) &i, sizeof( int ) );
          //outfile.write( (char*) &j, sizeof( int ) );
          //outfile.write( (char*) &k, sizeof( int ) );
          //outfile.write( (char*) &l, sizeof( int ) );
          outfile.write( (char*) &mo_index, sizeof( int ) );
          outfile.write( (char*) &val, sizeof( double ) );
 
          create_new_kl = false;

          create_first_cont = false;
          create_second_cont = false;
          create_third_cont = false;
          create_fourth_cont = false;
        }
      }
      }
    }
    }
    outfile.close(); 
}

double quarter_gamma_transform_new(
    const MatrixXd& evecsXd, const Eigen::MatrixXd& kernel_matr, const int& imo, const int& jmo, const int& kmo, const int& lmo,
    const bool& imo_is_new, const bool& jmo_is_new, const bool& kmo_is_new, const bool& lmo_is_new
){
    int nmo_scell = evecsXd.cols();
    static Eigen::VectorXd gamma_i     = Eigen::VectorXd::Zero(nmo_scell);
    static Eigen::VectorXd gamma_j     = Eigen::VectorXd::Zero(nmo_scell);
    static Eigen::VectorXd gamma_k     = Eigen::VectorXd::Zero(nmo_scell);
    static Eigen::VectorXd gamma_l     = Eigen::VectorXd::Zero(nmo_scell);
    static Eigen::VectorXd gamma_contr = Eigen::VectorXd::Zero(nmo_scell);
    double val = 0.0;
    if( imo_is_new || jmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_j( i ) = evecsXd( i, imo ) * evecsXd( i, jmo );
      //cout << "i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_j << endl; 
      gamma_contr = kernel_matr * gamma_j;
      //cout << "KERNEL i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_contr << endl; 
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd( i, kmo );
    }
    if( kmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd( i, kmo );
    }
    for( int i = 0; i < nmo_scell; ++i ) val += gamma_k( i ) * evecsXd( i, lmo );
    //if ( imo == 8 && kmo == 8 && jmo == 0 ) {
    //  cout << gamma_contr << endl;
    //  cout << " " << endl;
    //  cout << gamma_k << endl;
    //  cout << " " << endl;
    //}
    return val;
}

void gamma_moint_to_vector( 
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    const Eigen::MatrixXd& kernel_matr,
    const Eigen::MatrixXd& evecsXd,
    const int nmo_scell,
    size_t& arr_size,
    int pstart, int pend,
    int qstart, int qend,
    int rstart, int rend,
    int sstart, int send,
    int itol
){
    int ij_index, kl_index;
    size_t mo_index;
    double val, val2;
    double dtol = pow( 10., 1. * itol );
    
    size_t in_arr_size = indx_arr.size();
    int counter = 0;

    Eigen::VectorXd evgamma_i, evgamma_j, evgamma_k, evgamma_l, cont1, cont2, cont3;
    evgamma_i.resize( nmo_scell );
    evgamma_j.resize( nmo_scell );
    evgamma_k.resize( nmo_scell );
    evgamma_l.resize( nmo_scell );
    cont1.resize( nmo_scell );
    cont2.resize( nmo_scell );
    cont3.resize( nmo_scell );

    bool create_first_cont, create_second_cont, create_third_cont, create_fourth_cont;
    for( int i = pstart; i < pend; ++i ){
    evgamma_i = evecsXd.col( i );
    create_first_cont = true;
    for( int j = qstart; j < qend; ++j ){
      evgamma_j = evecsXd.col( j );
      create_second_cont = true;
      for( int k = rstart; k < rend; ++k ){
      evgamma_k = evecsXd.col( k );
      create_third_cont = true;
      for( int l = sstart; l < send; ++l ){
          evgamma_l = evecsXd.col( l );
          create_fourth_cont = true;
          mo_index = get_ijkl_index( i, j, k, l ); 
          //val = quarter_gamma_transform( i, j, k, l, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
          val = quarter_transform__gamma_point( evgamma_i, evgamma_j, evgamma_k, evgamma_l, cont1, cont2, cont3, kernel_matr, \
                                                nmo_scell, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
          //val = quarter_gamma_transform_new( evecsXd, kernel_matr, i, j, k, l, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont ); 
	  //printf ("%3d %3d %3d %3d   |   %20.16f  \n", i, j, k, l, val );
          //cout << "CREATING " << i << "," << j << "," << k << "," << l << " VAL : " << val << endl;

          if( fabs( val ) > dtol ){
            if( counter < in_arr_size ){
              dble_arr[ counter ] = val; 
              indx_arr[ counter ] = mo_index; 
            }else{
              dble_arr.push_back( val ); 
              indx_arr.push_back( mo_index ); 
            }
            counter++;
          }
          
          create_first_cont = false;
          create_second_cont = false;
          create_third_cont = false;
          create_fourth_cont = false;
      }
      }
    }
    }
    arr_size = counter; 
}

void moIntegralFactory::gamma_moint_to_vector( 
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    size_t& arr_size,
    int pstart, int pend,
    int qstart, int qend,
    int rstart, int rend,
    int sstart, int send,
    int itol
){
    int ij_index, kl_index;
    size_t mo_index;
    double val;
    double dtol = pow( 10., 1. * itol );
    
    size_t in_arr_size = indx_arr.size();
    int counter = 0;

    bool create_first_cont, create_second_cont, create_third_cont, create_fourth_cont;
    for( int i = pstart; i < pend; ++i ){
    create_first_cont = true;
    for( int j = qstart; j < qend; ++j ){
      create_second_cont = true;
      for( int k = rstart; k < rend; ++k ){
      create_third_cont = true;
        for( int l = sstart; l < send; ++l ){
        create_fourth_cont = true;
          mo_index = get_ijkl_index( i, j, k, l ); 
          val = quarter_gamma_transform( i, j, k, l, create_first_cont, create_second_cont, create_third_cont, create_fourth_cont );
          //cout << "CREATING " << mo_index << " VAL : " << val << endl;

          if( fabs( val ) > dtol ){
            if( counter < in_arr_size ){
              dble_arr[ counter ] = val; 
              indx_arr[ counter ] = mo_index; 
            }else{
              dble_arr.push_back( val ); 
              indx_arr.push_back( mo_index ); 
            }
            counter++;
          }
          
          create_first_cont = false;
          create_second_cont = false;
          create_third_cont = false;
          create_fourth_cont = false;
        }
      }
    }
    }
    arr_size = counter; 
}


double gamma_one_body_integrals(
    const Eigen::VectorXd& gamma_i, 
    const Eigen::VectorXd& gamma_j, 
    Eigen::VectorXd& temp_gamma_vec,
    const Eigen::MatrixXd& one_body_matr,
    int nmo_scell, 
    bool new_i 
){
    double val = 0.0;
    if( new_i ) temp_gamma_vec = one_body_matr * gamma_i;
    for( int i = 0; i < nmo_scell; ++i ) val += gamma_j( i ) * temp_gamma_vec( i );
    return val;
}


double gamma_fock_integral_ij(
    int p, int q,
    const Eigen::VectorXd& gamma_p, 
    const Eigen::VectorXd& gamma_q, 
    Eigen::VectorXd& gamma_contr_one_body,
    Eigen::VectorXd& gamma_contr1,
    Eigen::VectorXd& gamma_contr2,
    Eigen::VectorXd& gamma_contr3,
    const Eigen::MatrixXd& one_body_matr,
    const Eigen::MatrixXd& two_body_matr,
    const Eigen::MatrixXd& coulomb_matr,
    const Eigen::MatrixXd& evecsXd,
    double& xc_contribution,
    double& coulomb_contribution,
    const int nmo_scell,
    bool new_p_not_used,
    bool new_q_not_used
){
    bool new_i = true;
    bool new_p = true;
    bool new_q = true;
    double one_body_contribution = 0.0;
    double two_body_contribution = 0.0;
    double one_body_correction   = 0.0;
    int nocc = (int)(nmo_scell/2);
    Eigen::VectorXd gamma_i;
    gamma_i.resize( nmo_scell );
    for( int i = 0; i < nmo_scell; ++i ) gamma_i( i ) = 0.0;

    //cout << "FOCK ELEMENT (" << p << ", " << q << ") " << endl;

    one_body_contribution = gamma_one_body_integrals( gamma_p, gamma_q, gamma_contr_one_body, one_body_matr,
                                                      nmo_scell, new_p );

    new_p = true;
    new_q = true;
    double jint, xcint;
    xc_contribution = 0.0;
    coulomb_contribution = 0.0;
    for( int i = 0; i < nocc; ++i ){
      new_i = true;
      gamma_i = evecsXd.col( i ); 
      jint  = quarter_transform__gamma_point( gamma_p, gamma_q, gamma_i, gamma_i,
                gamma_contr1, gamma_contr2, gamma_contr3, two_body_matr, nmo_scell,
                new_p, new_q, new_i, new_i );
      xcint = quarter_transform__gamma_point( gamma_p, gamma_i, gamma_i, gamma_q,
                gamma_contr1, gamma_contr2, gamma_contr3, two_body_matr, nmo_scell,
                new_p, new_i, new_q, new_i );
      xc_contribution -= xcint;
      coulomb_contribution += 2. * jint;
      two_body_contribution += 2. * jint - xcint;
      
      //printf( "+2 CONTRIBUTION (%3d,%3d,%3d,%3d) : %20.16f \n", (p+1), (q+1), (i+1), (i+1), jint );
      //printf( "-1 CONTRIBUTION (%3d,%3d,%3d,%3d) : %20.16f \n", (p+1), (i+1), (q+1), (i+1), xcint );
    }


    //
    //
    // A separate loop makes sure the quarter transform behaves properly when the kernel is replaced
    //
    //
    new_p = true;
    new_q = true;
    for( int i = 0; i < nocc; ++i ){
      new_i = true;
      gamma_i = evecsXd.col( i ); 
      jint = 1. * quarter_transform__gamma_point( gamma_p, gamma_q, gamma_i, gamma_i, 
               gamma_contr1, gamma_contr2, gamma_contr3, coulomb_matr, nmo_scell,
               new_p, new_q, new_i, new_i );
      xcint = 1. * quarter_transform__gamma_point( gamma_p, gamma_q, gamma_i, gamma_i, 
               gamma_contr1, gamma_contr2, gamma_contr3, two_body_matr, nmo_scell,
               new_p, new_q, new_i, new_i );
      //printf( "CORRECTION (%3d,%3d,%3d,%3d) : %20.16f \n", (p+1), (q+1), (i+1), (i+1), jint );
      //printf( "CORRECTION (%3d,%3d,%3d,%3d) : %20.16f \n", (p+1), (q+1), (i+1), (i+1), xcint );
      one_body_correction += jint - xcint;
    }
    one_body_correction *= 2.;
    coulomb_contribution += one_body_correction; 

    double outval = one_body_contribution + two_body_contribution + one_body_correction;
    //printf( "FOCK ONE-BODY CONTRIBUTION %20.16f \n", one_body_contribution );
    //printf( "FOCK TWO-BODY CONTRIBUTION %20.16f \n", two_body_contribution );
    //printf( "FOCK ONE-BODY CORRECTION   %20.16f \n", one_body_correction   );
    //printf( "FOCK ELEMENT (%3d, %3d )      %20.16f \n", p, q, outval );
    return outval;
}

void check_gamma_fock(
    const Eigen::MatrixXd& one_body_matr,
    const Eigen::MatrixXd& two_body_matr,
    const Eigen::MatrixXd& coulomb_matr,
    const Eigen::MatrixXd& evecsXd,
    const Eigen::MatrixXd& evalsXd,
    const int nmo_scell
){
    Eigen::VectorXd gamma_contr_one_body, gamma_contr1, gamma_contr2, gamma_contr3, gamma_p, gamma_q; 
    gamma_contr_one_body.resize( nmo_scell );
    gamma_contr1.resize( nmo_scell );
    gamma_contr2.resize( nmo_scell );
    gamma_contr3.resize( nmo_scell );
    gamma_p.resize( nmo_scell );
    gamma_q.resize( nmo_scell );
    double fpq;
    bool new_p, new_q;
    double total = 0.0;
    int nocc = (int)(nmo_scell/2);
    double max_error = 0.0;
    double error = 0.0;
    double xc_contr, coulomb_contr;
    double exchange_energy, coulomb_energy;

    exchange_energy = 0.0;
    coulomb_energy = 0.0;
    
    for( int p = 0; p < nmo_scell; ++p ){
      new_p = true;
      gamma_p = evecsXd.col( p );
      for( int q = 0; q < nmo_scell; ++q ){
        new_q = true;
        gamma_q = evecsXd.col( q );
        //cout << "P " << gamma_p << endl;
        //cout << "Q " << gamma_q << endl;
        //cout << "KERNEL " << two_body_matr << endl;
        fpq = gamma_fock_integral_ij( p, q, gamma_p, gamma_q, gamma_contr_one_body, gamma_contr1, gamma_contr2, \
                                      gamma_contr3, one_body_matr, two_body_matr, coulomb_matr, evecsXd, xc_contr, coulomb_contr, nmo_scell , \
                                      new_p, new_q );
        if( p == q && p < (int)(nmo_scell/2) ){
          exchange_energy += xc_contr;
          coulomb_energy += coulomb_contr; 
        }
        if( p != q ) error = fabs( fpq - 0.0 ); 
        else error = fabs( fpq - evalsXd( p, 0 ) );

        if( error > max_error ) max_error = error;

        if( max_error > 1e-10 ){
	  printf( "ERROR : check_gamma_fock returned an error in the fock integrals of %20.16f \n", max_error );
	  printf( "        for fpq (p=%3d, q=%3d ).  Most likely this is not the lowest HF solution. \n", p, q );
	  printf( "        And with this, Brillouin's Theorem will not be satisfied.                 \n" );
          exit( EXIT_FAILURE );
	} 
        new_q = false;
        new_p = false;
      }
    }
    printf( "ERROR IN MO FOCK INTEGRALS / BRILLOUINS THEOREM = %20.14e \n", max_error );
    printf( "EXCHANGE ENERGY FROM MO's = %20.14e \n", exchange_energy );
    printf( "COULOMB  ENERGY FROM MO's = %20.14e \n", coulomb_energy );
}


void moIntegralFactory::Init(
    UnitCell& UCell,
    SuperCell& SCell,
    aoIntegralFactory& aoints,
    const PBC_CLASS& pbc,
    VMatrixXd& inevals,
    VMatrixXd& inevecsXd,
    VMatrixXcd& inevecsXcd
){
    printf( "SETTING UP INTEGRAL DRIVER... \n" );
    evals = inevals;

    nkpt = SCell.nkpt;
    ntrans = SCell.total_number_of_cells;

    nmo_ucell = UCell.nao;
    nmo_scell = SCell.nao;

    kernel_matr.resize( nmo_scell, nmo_scell );
    one_body_matr = Eigen::MatrixXd::Zero( nmo_scell, nmo_scell );
    coulomb_matr = Eigen::MatrixXd::Zero( nmo_scell, nmo_scell );

    for( int iao = 0; iao < nmo_ucell; ++iao ){
    for( int jao = 0; jao < nmo_ucell; ++jao ){
      for( int itrans = 0; itrans < ntrans; ++itrans ){
      for( int jtrans = 0; jtrans < ntrans; ++jtrans ){
        Vector3i tvec = SCell.reduced_t[ jtrans ] - SCell.reduced_t[ itrans ];
        bool found = false;
        int index = SCell.get_supercell_index( tvec( 0 ), tvec( 1 ), tvec( 2 ), &found );
        //
        //
        // making the one body matrix
        //
        //
        if( iao == jao && itrans == jtrans )
          one_body_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getOnsiteE( UCell.ElementStr[ iao ] );
        one_body_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getHopping( iao, jao, index );
        //
        //
        // Now adding a messy nuclear attraction contribution ....
        //
        //
        if( iao == jao && itrans == jtrans ){
 
        for( int jat = 0; jat < nmo_ucell; ++jat ){
          for( int ltrans = 0; ltrans < ntrans; ++ltrans ){
            Vector3i tvec1 = SCell.reduced_t[ ltrans ] - SCell.reduced_t[ itrans ];
            bool found1 = false;
            int index1 = SCell.get_supercell_index( tvec1( 0 ), tvec1( 1 ), tvec1( 2 ), &found1 );
              one_body_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) -= \
                         UCell.charge[ jat ] * aoints.getCoulombInt( iao, jat, index1 );
          }
        }

        } // end PPP condition that iao == jao
        //
        //
        // making the two body matrix
        //
        //
        kernel_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) = aoints.getXCInt( iao, jao, index );
        coulomb_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) = aoints.getCoulombInt( iao, jao, index );
        if( iao == jao && itrans == jtrans ){
          kernel_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getHubbardU( SCell.ElementStr[ iao ] );
          coulomb_matr( iao + itrans * nmo_ucell, jao + jtrans * nmo_ucell ) += aoints.getHubbardU( SCell.ElementStr[ iao ] );
        }
      }
      }
    }
    }

    Eigen::MatrixXd gamma_evecsXd = Eigen::MatrixXd::Zero( nmo_scell, nmo_scell );
    Eigen::MatrixXd gamma_evalsXd;
    gamma_evecsXd = inevecsXd.irrep( 0 );
    gamma_evalsXd = inevals.irrep( 0 );
    //check_gamma_fock( one_body_matr, kernel_matr, coulomb_matr, gamma_evecsXd, gamma_evalsXd, nmo_scell );

    //Eigen::VectorXd one_body_rhs, one_body_lhs, one_body_tvec;
    //one_body_rhs.resize( nmo_scell );
    //one_body_lhs.resize( nmo_scell );
    //one_body_tvec.resize( nmo_scell );
    //double one_body_energy = 0.0;
    //for( int i = 0; i < (int)(nmo_scell/2); ++i ){
    //  one_body_rhs = inevecsXd.irrep( 0 ).col( i );
    //  one_body_lhs = inevecsXd.irrep( 0 ).col( i );
    //  one_body_energy += gamma_one_body_integrals( one_body_lhs, one_body_rhs, one_body_tvec, one_body_matr, nmo_scell, true ); 
    //}
    //one_body_energy *= 2.;
    //printf( "ONE BODY ENERGY FROM MO INTEGRALS : %20.16f \n", ( one_body_energy / nkpt ) );

    if( pbc.type == GAMMA ){
      evecsXd = inevecsXd;
      gamma_contr.resize( nmo_scell );
      gamma_kl.resize( nmo_scell );
      gamma_ij.resize( nmo_scell );

      gamma_i.resize( nmo_scell );
      gamma_j.resize( nmo_scell );
      gamma_k.resize( nmo_scell );
      gamma_l.resize( nmo_scell );
    }else{
      printf( " - initializing cos and sin arguments... \n");
      evecsXcd = inevecsXcd;
      //setup_cos_sin( SCell ); 
    }

    double val_coulomb, val_xc, val_diff;
    val_coulomb = make_2eint_pqrs( coulomb_matr, evecsXd.irrep( 0 ), 0, 0, 0, 0 );
    printf( "VAL2 COULOMB : %20.16f \n", val_coulomb );
    val_xc = make_2eint_pqrs( kernel_matr, evecsXd.irrep( 0 ), 0, 0, 0, 0 );
    printf( "VAL2 XC      : %20.16f \n", val_xc );
    val_diff = fabs( val_coulomb - val_xc );
    printf( "DIFF         : %20.16f \n", val_diff );

    printf( "SETUP COMPLETE! \n" );
}

void moIntegralFactory::check_gamma_exchange_energy( std::vector< double >& gamma_ints, const double& ref_energy ){
    int nocc = (int)( nmo_scell / 2 );
    double exc = 0.0;
    for( int i = 0; i < nocc; ++i ){
    for( int j = 0; j < nocc; ++j ){ 
      //printf( "ADDING CONTRIBUTION FROM %3d   : %20.16f \n", get_ijkl_index( i, j, j, i ), gamma_ints[ get_ijkl_index( i, j, j, i ) ] );
      exc -= 1.0 * gamma_ints[ get_ijkl_index( i, j, i, j ) ];
    }
    }

    printf( "CHECKING EXCHANGE ENERGY MOINTEGRALS... \n" );
    printf( "  - EXCHANGE ENERGY WITH SCF         : %20.16f \n", ref_energy );
    printf( "  - EXCHANGE ENERGY WITH MO's        : %20.16f \n", exc );
    printf( "  - EXCHANGE ENERGY WITH MO's / nkpt : %20.16f \n", ( exc / nkpt ) );
}

double moIntegralFactory::quarter_gamma_transform(
    const int& imo, const int& jmo, const int& kmo, const int& lmo,
    const bool& imo_is_new, const bool& jmo_is_new, const bool& kmo_is_new, const bool& lmo_is_new
){
    double val = 0.0;
    if( imo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_j( i ) = evecsXd.irrep( 0 )( i, imo ) * evecsXd.irrep( 0 )( i, jmo );
      //cout << "i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_j << endl; 
      gamma_contr = kernel_matr * gamma_j;
      //cout << "KERNEL i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_contr << endl; 
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    if( jmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_j( i ) = evecsXd.irrep( 0 )( i, imo ) * evecsXd.irrep( 0 )( i, jmo );
      gamma_contr = kernel_matr * gamma_j;
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    if( kmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    for( int i = 0; i < nmo_scell; ++i ) val += gamma_k( i ) * evecsXd.irrep( 0 )( i, lmo );
    return val;
}

double moIntegralFactory::quarter_gamma_transform_coulomb(
    const int& imo, const int& jmo, const int& kmo, const int& lmo,
    const bool& imo_is_new, const bool& jmo_is_new, const bool& kmo_is_new, const bool& lmo_is_new
){
    double val = 0.0;
    if( imo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_j( i ) = evecsXd.irrep( 0 )( i, imo ) * evecsXd.irrep( 0 )( i, jmo );
      //cout << "i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_j << endl; 
      gamma_contr = coulomb_matr * gamma_j;
      //cout << "KERNEL i dot j " << imo << ", " << jmo << endl;
      //cout << gamma_contr << endl; 
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    if( jmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_j( i ) = evecsXd.irrep( 0 )( i, imo ) * evecsXd.irrep( 0 )( i, jmo );
      gamma_contr = coulomb_matr * gamma_j;
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    if( kmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_k( i ) = gamma_contr( i ) * evecsXd.irrep( 0 )( i, kmo );
    }
    for( int i = 0; i < nmo_scell; ++i ) val += gamma_k( i ) * evecsXd.irrep( 0 )( i, lmo );
    return val;
}

double quarter_transform__gamma_point(
    const Eigen::VectorXd& gamma_i, 
    const Eigen::VectorXd& gamma_j, 
    const Eigen::VectorXd& gamma_k, 
    const Eigen::VectorXd& gamma_l, 
    Eigen::VectorXd& gamma_contr1,
    Eigen::VectorXd& gamma_contr2,
    Eigen::VectorXd& gamma_contr3,
    const Eigen::MatrixXd& kernel_matr, 
    const int nmo_scell,
    const bool& imo_is_new, 
    const bool& jmo_is_new, 
    const bool& kmo_is_new, 
    const bool& lmo_is_new
){
    double val = 0.0;
    if( imo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_contr1( i ) = gamma_i( i ) * gamma_j( i );
      gamma_contr3 = kernel_matr * gamma_contr1;
      for( int i = 0; i < nmo_scell; ++i ) gamma_contr2( i ) = gamma_contr3( i ) * gamma_k( i ); 
    }
    if( jmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_contr1( i ) = gamma_i( i ) * gamma_j( i );
      gamma_contr3 = kernel_matr * gamma_contr1;
      for( int i = 0; i < nmo_scell; ++i ) gamma_contr2( i ) = gamma_contr3( i ) * gamma_k( i );
    }
    if( kmo_is_new ){
      for( int i = 0; i < nmo_scell; ++i ) gamma_contr2( i ) = gamma_contr3( i ) * gamma_k( i );
    }
    for( int i = 0; i < nmo_scell; ++i ) val += gamma_contr2( i ) * gamma_l( i );
    return val;
}

