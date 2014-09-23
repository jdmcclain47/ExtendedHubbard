#include <algorithm>
#include <vector>
#include <cstdio>
#include <complex>
#include <math.h>
#include <iostream>
#include <iomanip>
#include "assert.h"
#include "ewald.h"
#include "find_bounding_sphere.h"
#include "generate_tvecs.h"
#include "Eigen/Dense" 
using namespace std;
using namespace Eigen;

/* TODO : Look at CUTOFF ERRORS IN THE EWALD SUMMATION FORMULAE FOR POINT CHARGE SYSTEMS by Perram to get better a priori
 * error estimates for the ewald.  Right now, the errors are ~10^(-5) for most of these ewald functions, but one can
 * use the "_converge" option to get more energies within a given error tolerance
 */

static int suppressoutput = 0;
static double alpha;
static bool isalphaset = false;
static int nextra = 0;
static int nextra_energy = 0;

Eigen::Matrix3d create_Gmat( Eigen::Matrix3d& Tmat ){
    Matrix3d Gmat;
    double twopi = 2. * M_PI;
    RowVector3d abcross = Tmat.row( 1 ).cross( Tmat.row( 2 ) );
    Gmat.row( 0 ) = abcross / ( Tmat.row( 0 ).dot( abcross ) );

    abcross = Tmat.row( 2 ).cross( Tmat.row( 0 ) );
    Gmat.row( 1 ) = abcross / ( Tmat.row( 1 ).dot( abcross ) );

    abcross = Tmat.row( 0 ).cross( Tmat.row( 1 ) );
    Gmat.row( 2 ) = abcross / ( Tmat.row( 2 ).dot( abcross ) );

    Gmat *= twopi;

    return Gmat;
}


double dot_prod3( vector<double>& v1, vector<double>& v2 ){
   double dp3 = 0.0;
   for( int i = 0; i < v1.size(); i++ ){
      dp3 += v1[i]*v2[i];
   }
   return dp3;
}

void setalpha( double alpha_in ){
   alpha = alpha_in;
   isalphaset = true;
}

static Eigen::Matrix3d current_Tvec = Eigen::Matrix3d::Zero( 3, 3 );

void minimize_alpha_cost( Matrix3d& Tvec, Matrix3d& Gvec, double& step_size, vector< double >& alist, vector< size_t >& clist, int& current_step ){
  int nx, ny, nz;
  size_t gcost, rcost, totalcost;
  int max_steps = 1000;
  double calpha;
  if( current_step == 0 ) calpha = max( step_size, alpha - 0.5 * max_steps * step_size ); 
  else calpha = alist[ current_step - 1 ] + step_size;

  bool decreasing = true;

  double asqm4 = calpha * calpha * 4.;
  double gcutoff = sqrt( 20 * asqm4 ); 
  findBoundingSphere( Gvec, gcutoff, nx, ny, nz ); 
  gcost = nx * ny * nz;
  double rcutoff = 4.5 / calpha;
  findBoundingSphere( Tvec, rcutoff, nx, ny, nz ); 
  rcost = nx * ny * nz;
  totalcost = rcost + gcost;
  if( totalcost > (size_t)-1 || totalcost < 0 ) totalcost = (size_t)-1;

  std::string str_dec;
  if( current_step < max_steps ){
    if( current_step > 0 ) decreasing = ( totalcost <= clist[ current_step - 1 ] );
    str_dec = ( decreasing ) ? "TRUE" : "FALSE";
    printf( "CURRENT COST : %24lu     CURRENT ALPHA : %20.16f     DECREASING? %s \n", totalcost, calpha, str_dec.c_str() );
    if( decreasing ){
      alist.push_back( calpha );
      clist.push_back( totalcost );
      current_step++;
      minimize_alpha_cost( Tvec, Gvec, step_size, alist, clist, current_step );
    }else{
      return;
    }
  }
}


void setalpha( Matrix3d& Tvec, int nparticles, double volume ){
   double maxdist = 0.0;
   bool same_tmatrix = true;
   for( int i = 0; i < 3; ++i ) 
   for( int j = 0; j < 3; ++j )
     if( fabs( Tvec( i, j ) - current_Tvec( i, j ) ) > 1e-15 )
       same_tmatrix = false;
   if( same_tmatrix && isalphaset ) return;

   current_Tvec = Tvec;
   for( int i = 0; i < 3; i++ ){
      if( Tvec.row( i ).dot( Tvec.row( i ) ) > maxdist )
         maxdist = sqrt( Tvec.row( i ).dot( Tvec.row( i ) ) );
   }
   maxdist /= 2.0;
   alpha = 3.0 / maxdist; 
   Matrix3d Gvec = create_Gmat( Tvec );
   int current_step = 0;
   vector< double > alist;
   vector< size_t > clist;
   double step_size = 0.0005;
   minimize_alpha_cost( Tvec, Gvec, step_size, alist, clist, current_step );
   std::vector< size_t >::iterator iter = min_element( clist.begin(), clist.end() );
   int element_pos = ( iter - clist.begin() );
   alpha = alist[ element_pos ];
   isalphaset = true;
}

double get_alpha( Matrix3d& Tvec ){
    double out_alpha;
    out_alpha = alpha;
    return out_alpha; 
}

static int nextra_en_real = 0;
static int nextra_en_recip = 0;

double energy_ewald_converge_FAST( 
   Matrix3d& Tvec, 
   Matrix3d& Gvec, 
   std::vector< Eigen::Vector3d >& coord,
   vector< double >& charge,
   int itol
){
   // enEwald   : ewald energy
   // enLRange	: long-ranged part of ewald energy (reciprocal space)
   // enSRange	: short-ranged part of ewald energy (direct space)
   // enSI	: self-interaction part of ewald energy
   // enBK	: background part of ewald energy
   //                - to be used if a uniform background charge needs to
   //                  be subtracted from the input density to obtain the
   //                  ewald energy
   // alpha     : alpha term in ewald summation

   double enEwald  = 0.0;
   double enSRange = 0.0;
   double enLRange;
   double enSI     = 0.0;
   double enBK     = 0.0;

   double dtol = pow( 10., 1. * itol );

   //TODO find an alpha that minimized maxG + maxT

   double denom;

   // ... now creating the short-ranged term
   Vector3d T1 = Tvec.row(0);
   Vector3d T2 = Tvec.row(1);
   Vector3d T3 = Tvec.row(2);
   double volume = fabs( T1.dot( T2.cross( T3 ) ) );

   setalpha( Tvec, coord.size(), volume );

   if( !suppressoutput ){
     cout << "ALPHA : " << alpha << endl;
   }

   int nPoints = coord.size();
   double transvecLength;
   Vector3d transvec;
   Vector3d coord1;
   Vector3d coord2;
   Vector3d coordDiff;
   Vector3d G1 = Gvec.row(0);
   Vector3d G2 = Gvec.row(1);
   Vector3d G3 = Gvec.row(2);
   Vector3d recipvec;
   double gdotg;

   vector< double > xyzvals;
   vector< double > xyzvals_sq;
   int nvals;
   
   double maxG = 0.0;
   vector<int> mxG (3);
   vector<int> old_mxG (3, 0 );
   int nx, ny, nz;

   double old_enLRange = 1e20;
   complex<double> structureFactor;
   double asqm4 = alpha * alpha * 4.0;

   double gcutoff = sqrt( 20 * asqm4 ); 
   double rdotg;
   findBoundingSphere( Gvec, gcutoff, nx, ny, nz ); 

   cout << "========= EWALD ENERGY =========" << endl;
   cout << " CONVERGENCE = 10^(" << setw( 3 ) << itol << ")" << endl; 
   cout << "================================" << endl;

   //
   //
   // Starting convergence for the reciprocal space part...
   //
   //

   cout << "STARTING GSPACE CONVERGENCE..." << endl;
   enLRange = 0.0;
   double prefacLRange = (4.0 * M_PI / 2.0) / volume;
   while( nextra_en_recip == 0 || fabs( old_enLRange - enLRange ) > dtol ){

     old_enLRange = enLRange;

     mxG[ 0 ] = nx + nextra_en_recip;
     mxG[ 1 ] = ny + nextra_en_recip;
     mxG[ 2 ] = nz + nextra_en_recip;

     cout << "GSPACE CONV STEP " << setw( 3 ) << nextra_en_recip << flush;
     //cout << " - GENERATING GVECS FOR [ " << setw( 3 ) << old_mxG[ 0 ] << " : " << setw( 3 ) << mxG[ 0 ] << " ] ... " << flush;
     generate_tvecs( nvals, Gvec, old_mxG[ 0 ], mxG[ 0 ], xyzvals, xyzvals_sq );
     //cout << "done!" << endl;
     for( int i = 0; i < nvals; ++i ){
       structureFactor *= 0.0;
       for( int iat = 0; iat < nPoints; ++iat ){
         rdotg = xyzvals[ 3 * i + 0 ] * coord[ iat ]( 0 ) + \
                 xyzvals[ 3 * i + 1 ] * coord[ iat ]( 1 ) + \
                 xyzvals[ 3 * i + 2 ] * coord[ iat ]( 2 );
         real( structureFactor ) += charge[ iat ] * cos( rdotg );
         imag( structureFactor ) += charge[ iat ] * sin( rdotg );
       }
       gdotg = xyzvals_sq[ 3 * i + 0 ] + xyzvals_sq[ 3 * i + 1 ] + xyzvals_sq[ 3 * i + 2 ];
       //cout << "NEW : " << gdotg << " " << real( structureFactor ) << endl;
       if( gdotg > 1e-15 ){
          enLRange += prefacLRange * exp( - ( gdotg / asqm4) ) / gdotg * real( structureFactor * conj( structureFactor ) );
       }
     }

     printf( " - GPSACE ENERGY : %20.16f       DIFF : %20.16e \n", enLRange, fabs( enLRange - old_enLRange ) );
     fflush( stdout );

     old_mxG[ 0 ] = mxG[ 0 ]; 
     old_mxG[ 1 ] = mxG[ 1 ]; 
     old_mxG[ 2 ] = mxG[ 2 ]; 

     nextra_en_recip++;
   } // end convergence
   nextra_en_recip--;

   //
   //
   // Starting convergence for the reciprocal space part...
   //
   //
   vector< int > mxT ( 3 );
   vector< int > old_mxT ( 3, 0 );
   double rcutoff = 4.5 / alpha;
   double old_enSRange = 1e20;
   findBoundingSphere( Tvec, rcutoff, nx, ny, nz ); 

   cout << "STARTING RSPACE CONVERGENCE..." << endl;
   enSRange = 0.0;
   while( nextra_en_real == 0 || fabs( old_enSRange - enSRange ) > dtol ){

     old_enSRange = enSRange;

     mxT[ 0 ] = nx + nextra_en_real;
     mxT[ 1 ] = ny + nextra_en_real;
     mxT[ 2 ] = nz + nextra_en_real;

     cout << "RSPACE CONV STEP " << setw( 3 ) << nextra_en_real << flush;
     generate_tvecs( nvals, Tvec, old_mxT[ 0 ], mxT[ 0 ], xyzvals, xyzvals_sq );

     for( int i = 0; i < nvals; ++i ){
       for( int iat = 0; iat < nPoints; iat++ ){
          coord1 = coord[ iat ];
          for( int jat = 0; jat < nPoints; jat++ ){
             coord2 = coord[ jat ];
             coordDiff = coord1 - coord2;
             coordDiff( 0 ) += xyzvals[ 3 * i + 0 ];
             coordDiff( 1 ) += xyzvals[ 3 * i + 1 ];
             coordDiff( 2 ) += xyzvals[ 3 * i + 2 ];
             transvecLength = sqrt( coordDiff(0) * coordDiff(0) + coordDiff(1) * coordDiff(1) + coordDiff(2) * coordDiff(2) );
             if( transvecLength > 1e-10 ){
                enSRange += 0.5 * charge[iat] * charge[jat] * (erfc(alpha * transvecLength)/transvecLength);
             }
          }
       }

     }

     printf( " - RSPACE ENERGY : %20.16f       DIFF : %20.16e \n", enSRange, fabs( enSRange - old_enSRange ) );
     fflush( stdout );

     old_mxT[ 0 ] = mxT[ 0 ];
     old_mxT[ 1 ] = mxT[ 1 ];
     old_mxT[ 2 ] = mxT[ 2 ];

     nextra_en_real++;
   } // end convergence
   nextra_en_real--;

   if( !suppressoutput )
     printf( "LR TERM : %20.16f   summed over [ %3d, %3d, %3d ] \n", enLRange, mxG[ 0 ], mxG[ 1 ], mxG[ 2 ] );

   if( !suppressoutput )
     printf( "SR TERM : %20.16f   summed over [ %3d, %3d, %3d ] \n", enSRange, mxT[ 0 ], mxT[ 1 ], mxT[ 2 ] );

   //creating neutralizing background charge
   double totcharge = 0.0;
   for ( int iat = 0; iat < nPoints; iat++ ){
      totcharge += charge[ iat ];
   }
   enBK = - 0.5 * M_PI * totcharge * totcharge / alpha / alpha / volume;
   if( !suppressoutput )
     printf( "BK TERM : %20.16f \n", enBK );
   
   //creating self-interaction term due to use of ewald summation
   for( int iat = 0; iat < nPoints; iat++ ){
      enSI += charge[iat]*charge[iat];
   }
   enSI *= - alpha/sqrt(M_PI);
   if( !suppressoutput )
     printf( "SE TERM : %20.16f \n", enSI );

   enEwald = (enSRange + enLRange + enSI) + enBK;
   if( !suppressoutput )
     printf( "ENERGY  : %20.16f \n", enEwald );
   return enEwald;
}



double energy_ewald( Matrix3d& Tvec, Matrix3d& Gvec, Matrix< double, Eigen::Dynamic, 3 >& coord, vector< double >& charge ){

   cout << fixed << setprecision( 16 );

   // enEwald   : ewald energy
   // enLRange	: long-ranged part of ewald energy (reciprocal space)
   // enSRange	: short-ranged part of ewald energy (direct space)
   // enSI	: self-interaction part of ewald energy
   // enBK	: background part of ewald energy
   //                - to be used if a uniform background charge needs to
   //                  be subtracted from the input density to obtain the
   //                  ewald energy
   // alpha     : alpha term in ewald summation

   double enEwald = 0.0;
   double enLRange = 0.0;
   double enSRange = 0.0;
   double enSI     = 0.0;
   double enBK     = 0.0;

   //TODO find an alpha that minimized maxG + maxT

   double denom;

   // ... now creating the short-ranged term
   Vector3d T1 = Tvec.row(0);
   Vector3d T2 = Tvec.row(1);
   Vector3d T3 = Tvec.row(2);
   double volume = fabs( T1.dot( T2.cross( T3 ) ) );

      setalpha( Tvec, coord.rows(), volume );

   if( !suppressoutput ){
     cout << "ALPHA : " << alpha << endl;
   }

   int nPoints = coord.rows();
   double transvecLength;
   Vector3d transvec;
   Vector3d coord1;
   Vector3d coord2;
   Vector3d coordDiff;
   Vector3d G1 = Gvec.row(0);
   Vector3d G2 = Gvec.row(1);
   Vector3d G3 = Gvec.row(2);
   Vector3d recipvec;
   double gdotg;
   
   // finding Gvectors to sum over for when argument in exponent is less than 1e-09 = exp(-20)
   // So we want to find where 
   //                    n^2 |G|^2 / alpha^2 / 4.0 > 20.0
   // or in other words, where
   //                    n > sqrt (20.0 * 4.0 * alpha^2 / |G|^2 )
   double maxG = 0.0;
   vector<int> mxG (3);
   mxG[ 0 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G1.dot( G1 ) ) ) + nextra_energy;
   mxG[ 1 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G2.dot( G2 ) ) ) + nextra_energy;
   mxG[ 2 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G3.dot( G3 ) ) ) + nextra_energy;

   complex<double> structureFactor;
   for(int ix=-mxG[0];ix<mxG[0]+1;ix++){
      for(int iy=-mxG[1];iy<mxG[1]+1;iy++){
         for(int iz=-mxG[2];iz<mxG[2]+1;iz++){
            for(int idim = 0;idim < 3;idim++){
               recipvec(idim) = 1.0*ix*G1(idim) + 1.0*iy*G2(idim) + 1.0*iz*G3(idim); 
            }
            structureFactor *= 0.0;
            for( int iat = 0; iat < nPoints; iat++ ){
               real(structureFactor) += charge[ iat ] * cos( recipvec.dot( coord.row( iat ) ) );
               imag(structureFactor) += charge[ iat ] * sin( recipvec.dot( coord.row( iat ) ) );
            } 
            //cout << structureFactor << endl;
            gdotg = recipvec.dot( recipvec );
            //cout << "OLD : " << gdotg << " " << real( structureFactor ) << endl;
            if( abs( ix ) + abs( iy ) + abs( iz ) != 0 ){
               enLRange += exp( -(gdotg/alpha/alpha/4.0) ) / gdotg * real( structureFactor * conj( structureFactor ) );
            }
         }
      }
   }
   // factor of 1/2 to account for double counting
   enLRange *= (4.0 * M_PI / 2.0) / volume;
   if( !suppressoutput )
   cout << "LR term " << enLRange << " summed over " << mxG[0] << " " << mxG[1] << " " << mxG[2] << endl;



   // finding tvectors to sum over for when argument in erfc is less than 4.5, giving 1e-10 tolerance 
   // In doing so, we have that the argument inside the erfc() is given by
   //               erfc( alpha * n * T)
   // where n is the number of translations in a given x,y,z direction.  Thus we want to have
   //               alpha * n * T < 4.5
   double maxT = 0.0;
   vector<int> mxT (3);
   mxT[ 0 ] = ceil ( 4.5 / alpha / sqrt( T1.dot( T1 ) ) ) + nextra_energy;
   mxT[ 1 ] = ceil ( 4.5 / alpha / sqrt( T2.dot( T2 ) ) ) + nextra_energy;
   mxT[ 2 ] = ceil ( 4.5 / alpha / sqrt( T3.dot( T3 ) ) ) + nextra_energy; 

   for(int ix=-mxT[0];ix<mxT[0]+1;ix++){
      for(int iy=-mxT[1];iy<mxT[1]+1;iy++){
         for(int iz=-mxT[2];iz<mxT[2]+1;iz++){

            for(int idim = 0;idim < 3;idim++){
               transvec(idim) = 1.0 * ix * T1(idim) + 1.0 * iy * T2(idim) + 1.0 * iz * T3(idim); 
            }

            for( int iat = 0; iat < nPoints; iat++ ){
               coord1 = coord.row( iat );
               for( int jat = 0; jat < nPoints; jat++ ){
                  coord2 = coord.row( jat );
                  coordDiff = coord1 - coord2 + transvec; 
                  transvecLength = sqrt( coordDiff(0) * coordDiff(0) + coordDiff(1) * coordDiff(1) + coordDiff(2) * coordDiff(2) );
                  if( transvecLength > 1e-10 ){
                     enSRange += charge[iat] * charge[jat] * (erfc(alpha * transvecLength)/transvecLength);
                  }
               }
            }

         }
      }
   }
   enSRange *= 0.5;

   cout << setprecision(16);
   if( !suppressoutput )
   cout << "SR term " << setw(16) << enSRange << " summed over " << mxT[0] << " " << mxT[1] << " " << mxT[2] << endl;

   //creating neutralizing background charge
   double totcharge = 0.0;
   for ( int iat = 0; iat < nPoints; iat++ ){
      totcharge += charge[ iat ];
   }
   enBK = - 0.5 * M_PI * totcharge * totcharge / alpha / alpha / volume;
   if( !suppressoutput )
   cout << "BK term " << setw(16) << enBK << endl;
   //creating self-interaction term due to use of ewald summation
   for( int iat = 0; iat < nPoints; iat++ ){
      enSI += charge[iat]*charge[iat];
   }
   enSI *= - alpha/sqrt(M_PI);
   if( !suppressoutput )
   cout << "SE term " << setw(16) << enSI << endl;

   enEwald = (enSRange + enLRange + enSI) + enBK;
   cout << setprecision( 16 );
   if( !suppressoutput )
   cout << enEwald << endl;
   return enEwald;
}

double energy_ewald_converge(
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec,
    std::vector< Eigen::Vector3d >& coord,
    std::vector<double>& charge,
    int itol
){
    double dtol = pow( 10., 1. * itol );
    Eigen::Matrix< double, Eigen::Dynamic, 3 > newcoord;
    newcoord.resize( coord.size(), 3 );
    for( int i = 0; i < coord.size(); ++i ){
        newcoord( i, 0 ) = coord[ i ]( 0 );
        newcoord( i, 1 ) = coord[ i ]( 1 );
        newcoord( i, 2 ) = coord[ i ]( 2 );
    }
    double ewald_energy = energy_ewald( Tvec, Gvec, newcoord, charge );
    double old_ewald_energy = ewald_energy + 999.;
    cout << fixed << setprecision( 16 );
    cout << "EWALD CONVERGING ENERGY (NUCLEAR)..." << endl;
    while( fabs( ewald_energy - old_ewald_energy ) > dtol ){
      old_ewald_energy = ewald_energy;
      nextra_energy++;
      ewald_energy =  energy_ewald( Tvec, Gvec, newcoord, charge );
      cout << " - EXTRA STEPS " << nextra_energy << " OLD ENERGY = " << old_ewald_energy << "   NEW ENERGY = " << ewald_energy;
      cout << "   DIFF = " << fabs( ewald_energy - old_ewald_energy ) << endl;
    }
    nextra_energy--;
    return ewald_energy;
}



double energy_ewald(
    Eigen::Matrix3d& Tvec,
    Eigen::Matrix3d& Gvec,
    std::vector< Eigen::Vector3d >& coord,
    std::vector<double>& charge
){
    Eigen::Matrix< double, Eigen::Dynamic, 3 > newcoord;
    newcoord.resize( coord.size(), 3 );
    for( int i = 0; i < coord.size(); ++i ){
        newcoord( i, 0 ) = coord[ i ]( 0 );
        newcoord( i, 1 ) = coord[ i ]( 1 );
        newcoord( i, 2 ) = coord[ i ]( 2 );
    }
    return energy_ewald( Tvec, Gvec, newcoord, charge );
}

double potential_ewald( Vector3d pos, Matrix3d& Tvec, Matrix3d& Gvec, std::vector< Eigen::Vector3d >& coord, vector< double >& charge ){
    Eigen::Matrix< double, Eigen::Dynamic, 3 > newcoord;
    newcoord.resize( coord.size(), 3 );
    for( int i = 0; i < coord.size(); ++i ){
        newcoord( i, 0 ) = coord[ i ]( 0 );
        newcoord( i, 1 ) = coord[ i ]( 1 );
        newcoord( i, 2 ) = coord[ i ]( 2 );
    }
    VectorXd pos3d = pos;
    return potential_ewald( pos3d, Tvec, Gvec, newcoord, charge ); 
}


double potential_ewald( Vector3d pos, Matrix3d& Tvec, Matrix3d& Gvec, Matrix< double, Eigen::Dynamic, 3 >& coord, vector< double >& charge ){
   double potEwald;
   VectorXd pos3d = pos;
   potEwald = potential_ewald( pos3d, Tvec, Gvec, coord, charge );
   return potEwald;
}

double potential_ewald_converge( Vector3d pos, Matrix3d& Tvec, Matrix3d& Gvec, Vector3d& coord, vector< double >& charge, int itol ){
   double dtol = pow( 10., 1. * itol );
   double old_ewald_energy = 999.9;
   double ewald_energy = potential_ewald( pos, Tvec, Gvec, coord, charge );
   while( fabs( old_ewald_energy - ewald_energy ) > dtol ){
     old_ewald_energy = ewald_energy;
     nextra++;
     ewald_energy = potential_ewald( pos, Tvec, Gvec, coord, charge );
   }
   nextra--;
   return ewald_energy;
}

double potential_ewald( Vector3d pos, Matrix3d& Tvec, Matrix3d& Gvec, Vector3d& coord, vector< double >& charge ){ 
   double potEwald;
   Matrix< double, Eigen::Dynamic, 3 > coordmat;
   coordmat.resize( 1, 3 );
   coordmat( 0, 0 ) = coord(0);
   coordmat( 0, 1 ) = coord(1);
   coordmat( 0, 2 ) = coord(2);
   potEwald = potential_ewald( pos, Tvec, Gvec, coordmat, charge );
   return potEwald;
}

double self_potential_ewald_converge( Matrix3d& Tvec, Matrix3d& Gvec, int itol ){
   Vector3d pos1, pos2;
   pos1 << 0.0, 0.0, 0.0; 
   pos2 << 0.0, 0.0, 0.0; 
   vector< double > charge ( 1, 1.0 );
   double dtol = pow( 10., 1. * itol );
   double old_ewald_energy = 999.9;
   double ewald_energy = potential_ewald( pos1, Tvec, Gvec, pos2, charge );
   cout << fixed << setprecision( 16 );
   cout << "EWALD CONVERGING SELF POTENTIAL..." << endl;
   while( fabs( old_ewald_energy - ewald_energy ) > dtol ){
     old_ewald_energy = ewald_energy;
     nextra++;
     ewald_energy = potential_ewald( pos1, Tvec, Gvec, pos2, charge );
     cout << " - EXTRA STEPS " << nextra << " OLD ENERGY = " << old_ewald_energy << "   NEW ENERGY = " << ewald_energy;
     cout << "   DIFF = " << fabs( ewald_energy - old_ewald_energy ) << endl;
   }
   nextra--;
   return ewald_energy; 
}

double self_potential_ewald( Matrix3d& Tvec, Matrix3d& Gvec ){
   Vector3d pos1, pos2;
   pos1 << 0.0, 0.0, 0.0; 
   pos2 << 0.0, 0.0, 0.0; 
   vector< double > charge ( 1, 1.0 );
   double potEwald = potential_ewald( pos1, Tvec, Gvec, pos2, charge );
   return potEwald;
}

static int nextra_pot_real = 0;
static int nextra_pot_recip = 0;

double self_potential_ewald_converge_FAST(
   Eigen::Matrix3d& Tvec, 
   Eigen::Matrix3d& Gvec, 
   int itol
){
   std::vector< double > charge (1, 1.0 );
   std::vector< Eigen::Vector3d > coord;
   Eigen::Vector3d zerovec;
   zerovec( 0 ) = 0.0;
   zerovec( 1 ) = 0.0;
   zerovec( 2 ) = 0.0;
   coord.push_back( zerovec );
   double val; 
   val = potential_ewald_converge_FAST( zerovec, Tvec, Gvec, coord, charge, itol );
   return val;
};

double potential_ewald_converge_FAST( 
   Eigen::Vector3d& pos, 
   Matrix3d& Tvec, 
   Matrix3d& Gvec, 
   Eigen::Vector3d& coord1, 
   std::vector< double >& charge,
   int itol
){
   vector< Eigen::Vector3d > coord;
   coord.push_back( coord1 );
   double val = potential_ewald_converge_FAST( pos, Tvec, Gvec, coord, charge, itol );
   return val;
}

double potential_ewald_converge_FAST( 
   Eigen::Vector3d& pos, 
   Matrix3d& Tvec, 
   Matrix3d& Gvec, 
   std::vector< Eigen::Vector3d >& coord, 
   std::vector< double >& charge,
   int itol
){
 

   double dtol = std::pow( 10., 1. * itol );

   // potEwald   : ewald energy
   // potLRange	: long-ranged part of ewald energy (reciprocal space)
   // potSRange	: short-ranged part of ewald energy (direct space)
   // potSI	: self-interaction part of ewald energy
   // potBK	: background part of ewald energy
   //                - to be used if a uniform background charge needs to
   //                  be subtracted from the input density to obtain the
   //                  ewald energy
   // alpha     : alpha term in ewald summation

   double potEwald  = 0.0;
   double potLRange = 0.0;
   double potSRange = 0.0;
   double potSI     = 0.0;
   double potBK     = 0.0;

   Vector3d T1 = Tvec.row(0);
   Vector3d T2 = Tvec.row(1);
   Vector3d T3 = Tvec.row(2);
   double volume = fabs( T1.dot(T2.cross(T3)) );
   //
   //
   // ... now creating the short-ranged term
   //
   //

      setalpha( Tvec, coord.size(), volume );

   int nvals;
   int nPoints = coord.size();
   double transvecLength;
   Vector3d transvec, coord1, coord2, coordDiff;
   Vector3d recipvec, recipvecx, recipvecy;
   double gdotg;
   
   double maxG = 0.0;
   vector< int > mxG (3);
   vector< int > old_mxG ( 3, 0 );
   int nx, ny, nz;

   double old_potLRange = 1e20;
   complex<double> structureFactor;
   double asqm4 = alpha * alpha * 4.0;

   double gcutoff = sqrt( 20. * asqm4 ); 
   double rdotg;
   findBoundingSphere( Gvec, gcutoff, nx, ny, nz ); 
   nx = max( nx, 1 );
   ny = max( ny, 1 );
   nz = max( nz, 1 );

   cout << "========= EWALD POTENTIAL =========" << endl;
   cout << " CONVERGENCE = 10^(" << setw( 3 ) << itol << ")" << endl; 
   cout << " ALPHA VALUE = " << setw( 20 ) << setprecision( 16 ) << alpha << endl; 
   cout << "===================================" << endl;

   std::vector< double > xyzvals;
   std::vector< double > xyzvals_sq;

   //
   //
   // Starting convergence for the reciprocal space part...
   //
   //

   cout << "STARTING GSPACE CONVERGENCE..." << endl;
   potLRange = 0.0;
   double prefacLRange = (4.0 * M_PI ) / volume;
   while( nextra_pot_recip == 0 || fabs( old_potLRange - potLRange ) > dtol ){

     old_potLRange = potLRange;

     mxG[ 0 ] = nx + nextra_pot_recip;
     mxG[ 1 ] = ny + nextra_pot_recip;
     mxG[ 2 ] = nz + nextra_pot_recip;

     cout << "GSPACE CONV STEP " << setw( 3 ) << nextra_pot_recip << flush;
     //cout << " - GENERATING GVECS FOR [ " << setw( 3 ) << old_mxG[ 0 ] << " : " << setw( 3 ) << mxG[ 0 ] << " ] ... " << flush;
     generate_tvecs( nvals, Gvec, old_mxG[ 0 ], mxG[ 0 ], xyzvals, xyzvals_sq );
     //cout << "done!" << endl;
     for( int i = 0; i < nvals; ++i ){
       structureFactor *= 0.0;
       for( int iat = 0; iat < nPoints; ++iat ){
         rdotg = xyzvals[ 3 * i + 0 ] * ( coord[ iat ]( 0 ) - pos( 0 ) ) + \
                 xyzvals[ 3 * i + 1 ] * ( coord[ iat ]( 1 ) - pos( 1 ) ) + \
                 xyzvals[ 3 * i + 2 ] * ( coord[ iat ]( 2 ) - pos( 2 ) );
         real( structureFactor ) += charge[ iat ] * cos( rdotg );
       }
       gdotg = xyzvals_sq[ 3 * i + 0 ] + xyzvals_sq[ 3 * i + 1 ] + xyzvals_sq[ 3 * i + 2 ];
       //cout << "NEW : " << gdotg << " " << real( structureFactor ) << endl;
       if( gdotg > 1e-15 ){
          potLRange += prefacLRange * exp( - ( gdotg / asqm4) ) / gdotg * real( structureFactor );
       }
     }

     printf( " - GPSACE ENERGY : %20.16f       DIFF : %20.16e \n", potLRange, fabs( potLRange - old_potLRange ) );
     fflush( stdout );

     old_mxG[ 0 ] = mxG[ 0 ]; 
     old_mxG[ 1 ] = mxG[ 1 ]; 
     old_mxG[ 2 ] = mxG[ 2 ]; 

     nextra_pot_recip++;
   } // end convergence
   nextra_pot_recip--;
   nextra_pot_recip--;
   if( nextra_pot_recip < 0 ) nextra_pot_recip = 0;

   //
   //
   // Starting convergence for the reciprocal space part...
   //
   //
   vector< int > mxT ( 3 );
   vector< int > old_mxT ( 3, 0 );
   double rcutoff = 4.5 / alpha;
   double old_potSRange = 1e20;
   findBoundingSphere( Tvec, rcutoff, nx, ny, nz ); 

   cout << "STARTING RSPACE CONVERGENCE..." << endl;
   potSRange = 0.0;
   while( nextra_pot_real == 0 || fabs( old_potSRange - potSRange ) > dtol ){

     old_potSRange = potSRange;

     mxT[ 0 ] = nx + nextra_pot_real;
     mxT[ 1 ] = ny + nextra_pot_real;
     mxT[ 2 ] = nz + nextra_pot_real;

     cout << "RSPACE CONV STEP " << setw( 3 ) << nextra_pot_real << flush;
     generate_tvecs( nvals, Tvec, old_mxT[ 0 ], mxT[ 0 ], xyzvals, xyzvals_sq );

     for( int i = 0; i < nvals; ++i ){
       for( int iat = 0; iat < nPoints; iat++ ){
          coordDiff( 0 ) = xyzvals[ 3 * i + 0 ] + coord[ iat ]( 0 ) - pos( 0 );
          coordDiff( 1 ) = xyzvals[ 3 * i + 1 ] + coord[ iat ]( 1 ) - pos( 1 );
          coordDiff( 2 ) = xyzvals[ 3 * i + 2 ] + coord[ iat ]( 2 ) - pos( 2 );
          transvecLength = sqrt( coordDiff(0) * coordDiff(0) + coordDiff(1) * coordDiff(1) + coordDiff(2) * coordDiff(2) );
          if( transvecLength > 1e-10 ){
             potSRange += charge[iat] * (erfc(alpha * transvecLength)/transvecLength);
          }
       }

     }

     printf( " - RSPACE ENERGY : %20.16f       DIFF : %20.16e \n", potSRange, fabs( potSRange - old_potSRange ) );
     fflush( stdout );

     old_mxT[ 0 ] = mxT[ 0 ];
     old_mxT[ 1 ] = mxT[ 1 ];
     old_mxT[ 2 ] = mxT[ 2 ];

     nextra_pot_real++;
   } // end convergence
   nextra_pot_real--;
   nextra_pot_real--;
   if( nextra_pot_real < 0 ) nextra_pot_real = 0;

   if( !suppressoutput )
     printf( "LR TERM : %20.16f   summed over [ %3d, %3d, %3d ] \n", potLRange, mxG[ 0 ], mxG[ 1 ], mxG[ 2 ] );

   if( !suppressoutput )
     printf( "SR TERM : %20.16f   summed over [ %3d, %3d, %3d ] \n", potSRange, mxT[ 0 ], mxT[ 1 ], mxT[ 2 ] );


   //creating neutralizing background charge
   double totcharge = 0.0;
   for ( int iat = 0; iat < nPoints; iat++ ){
      totcharge += charge[ iat ];
   }
   potBK = - 1.0 * M_PI * totcharge * totcharge / alpha / alpha / volume;
   if( !suppressoutput )
     printf( "BK TERM : %20.16f \n", potBK );
   
   // see if any atom has that same coordinate for self interaction term
   potSI = 0.0;
   for( int iat = 0; iat < nPoints; iat++ ){
      coord1(0) = coord[ iat ]( 0 ) - pos(0);
      coord1(1) = coord[ iat ]( 1 ) - pos(1);
      coord1(2) = coord[ iat ]( 2 ) - pos(2);
      if ( coord1.dot( coord1 ) < 1e-10 ){
         potSI += charge[iat];
      }
   }
   potSI *= - ( 2.0 ) * alpha/sqrt(M_PI);
   if( !suppressoutput )
     printf( "SE TERM : %20.16f \n", potSI );

   potEwald = (potSRange + potLRange + potSI) + potBK;
   if( !suppressoutput )
     printf( "POTEN.  : %20.16f \n", potEwald );
   return potEwald;
}

double potential_ewald( 
   VectorXd pos, 
   Matrix3d& Tvec, 
   Matrix3d& Gvec, 
   Matrix< double, Eigen::Dynamic, 3 >& coord, 
   vector< double >& charge 
){
 
   if( !suppressoutput ){ 
      cout << "BEGINNING EWALD" << endl;
      cout << fixed << setprecision(16);
      cout << "ALPHA : " << alpha << endl;
   }


   // potEwald   : ewald energy
   // potLRange	: long-ranged part of ewald energy (reciprocal space)
   // potSRange	: short-ranged part of ewald energy (direct space)
   // potSI	: self-interaction part of ewald energy
   // potBK	: background part of ewald energy
   //                - to be used if a uniform background charge needs to
   //                  be subtracted from the input density to obtain the
   //                  ewald energy
   // alpha     : alpha term in ewald summation

   double potEwald = 0.0;
   double potLRange = 0.0;
   double potSRange = 0.0;
   double potSI     = 0.0;
   double potBK     = 0.0;

   //TODO find an alpha that minimized maxG + maxT

   double denom;

   // ... now creating the short-ranged term
   Vector3d T1 = Tvec.row(0);
   Vector3d T2 = Tvec.row(1);
   Vector3d T3 = Tvec.row(2);
   double volume = fabs( T1.dot(T2.cross(T3)) );

      setalpha( Tvec, coord.rows(), volume );

   int nPoints = coord.rows();
   double transvecLength;
   Vector3d transvec;
   Vector3d coord1;
   Vector3d coord2;
   Vector3d coordDiff;
   //finding Gvectors to sum over for when argument in exponent is less than 1e-09 = exp(-20)
   Vector3d G1 = Gvec.row(0);
   Vector3d G2 = Gvec.row(1);
   Vector3d G3 = Gvec.row(2);
   Vector3d recipvec, recipvecx, recipvecy;
   double gdotg;
   
   double maxG = 0.0;
   vector<int> mxG (3);
   mxG[ 0 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G1.dot( G1 ) ) ) + nextra;
   mxG[ 1 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G2.dot( G2 ) ) ) + nextra;
   mxG[ 2 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G3.dot( G3 ) ) ) + nextra;

   double structureFactor;
   Vector3d relcoord;
   double exparg;
   for(int ix=-mxG[0];ix<mxG[0]+1;ix++){
      recipvecx = ix * G1;
      for(int iy=-mxG[1];iy<mxG[1]+1;iy++){
         recipvecy = iy * G2;
         for(int iz=-mxG[2];iz<mxG[2]+1;iz++){
            recipvec = recipvecx + recipvecy + iz * G3;
            gdotg = recipvec.dot( recipvec );
            exparg = - gdotg/alpha/alpha/4.0;
            structureFactor *= 0.0;
            if( abs( ix ) + abs( iy ) + abs( iz ) != 0 ){
               for( int iat = 0; iat < nPoints; iat++ ){
                  relcoord(0) = coord( iat, 0 ) - pos(0);
                  relcoord(1) = coord( iat, 1 ) - pos(1);
                  relcoord(2) = coord( iat, 2 ) - pos(2);
                  structureFactor += charge[ iat ] * cos( recipvec.dot( relcoord ) ) * exp( exparg )/gdotg ;
                  //cout << structureFactor << " made from " << cos( recipvec.dot( relcoord ) ) << " " << exp( -(gdotg/alpha/alpha/4.0) ) << " " << gdotg << endl;
               }
            } 
            potLRange +=  structureFactor ;
         }
      }
   }
   potLRange *= ( 4.0 ) * M_PI / volume;

   if( !suppressoutput ){
      cout << "LR term " << potLRange << " summed over " << mxG[0] << " " << mxG[1] << " " << mxG[2] << endl;
   }


   // finding tvectors to sum over for when argument in erfc is less than 4.5, giving 1e-10 tolerance 
   // In doing so, we have that the argument inside the erfc() is given by
   //               erfc( alpha * n * T)
   // where n is the number of translations in a given x,y,z direction.  Thus we want to have
   //               alpha * n * T < 4.5
   double maxT = 0.0;
   vector<int> mxT (3);
   mxT[ 0 ] = ceil ( 4.5 / alpha / sqrt( T1.dot( T1 ) ) ) + nextra;
   mxT[ 1 ] = ceil ( 4.5 / alpha / sqrt( T2.dot( T2 ) ) ) + nextra;
   mxT[ 2 ] = ceil ( 4.5 / alpha / sqrt( T3.dot( T3 ) ) ) + nextra; 

   for(int ix=-mxT[0];ix<mxT[0]+1;ix++){
      for(int iy=-mxT[1];iy<mxT[1]+1;iy++){
         for(int iz=-mxT[2];iz<mxT[2]+1;iz++){

            for(int idim = 0;idim < 3;idim++){
               transvec(idim) = 1.0 * ix * T1(idim) + 1.0 * iy * T2(idim) + 1.0 * iz * T3(idim); 
            }

            for( int iat = 0; iat < nPoints; iat++ ){
               coord1(0) = coord( iat, 0 ) - pos(0);
               coord1(1) = coord( iat, 1 ) - pos(1);
               coord1(2) = coord( iat, 2 ) - pos(2);
               coordDiff = coord1 + transvec; 
               transvecLength = sqrt( coordDiff.dot( coordDiff ) );
               if( transvecLength > 1e-10 ){
                  potSRange += charge[iat] * (erfc(alpha * transvecLength)/transvecLength);
               }
            }

         }
      }
   }

   //cout << setprecision(8);
   if( !suppressoutput ){
      cout << "SR term " << setw(12) << potSRange << " summed over " << mxT[0] << " " << mxT[1] << " " << mxT[2] << endl;
   }

   //creating neutralizing background charge
   double totcharge = 0.0;
   for ( int iat = 0; iat < nPoints; iat++ ){
      totcharge += charge[ iat ];
   }
   potBK = - 1.0 * M_PI * totcharge / alpha / alpha / volume;
   //potBK = 0.0; 

   if( !suppressoutput ){
      cout << "BK term " << setw(12) << potBK << endl;
   }
   // seeing whether we are calculating the potential where an ion is located, in which case we need to subtract out
   // a self interaction term from the long ranged potential

   //    first we translate to reference cell
   /* THE FOLLOWING ONLY WORKS FOR ORTHOGONAL LATTICE VECTORS... */
   /*
   for( int i = 0; i < 3; i++ ){
      double tempdotprod = 0.0;
      for( int j = 0; j < 3; j++ ){
         tempdotprod += pos( j ) * Tvec( i, j );
      }
      tempdotprod /= Tvec.row( i ).dot( Tvec.row( i ) );
      if( tempdotprod > 1.0 - 1e-04 ){
         pos( 0 ) -= int( tempdotprod ) * Tvec( i, 0 );
         pos( 1 ) -= int( tempdotprod ) * Tvec( i, 1 );
         pos( 2 ) -= int( tempdotprod ) * Tvec( i, 2 );
      } 
      if( tempdotprod < 0.0 + 1e-04 ){
         pos( 0 ) += int( tempdotprod ) * Tvec( i, 0 );
         pos( 1 ) += int( tempdotprod ) * Tvec( i, 1 );
         pos( 2 ) += int( tempdotprod ) * Tvec( i, 2 );
      } 
   }
   */

   //    see if any atom has that same coordinate 
   for( int iat = 0; iat < nPoints; iat++ ){
      coord1(0) = coord( iat, 0 ) - pos(0);
      coord1(1) = coord( iat, 1 ) - pos(1);
      coord1(2) = coord( iat, 2 ) - pos(2);
      if ( coord1.dot( coord1 ) < 1e-10 ){
         potSI += charge[iat];
      }
   }
   potSI *= - ( 2.0 ) * alpha/sqrt(M_PI);
   if( !suppressoutput ){
      cout << "SE term " << setw(12) << potSI << endl;
   }

   potEwald =  (potSRange + potLRange + potSI) + potBK;
//   potEwald =  (potSRange + potLRange + potSI);
   if( !suppressoutput ){
      cout << potEwald << endl;
   }
   return potEwald;
}


double self_potential_ewald__2D(
    Matrix3d& Tvec,
    Matrix3d& Gvec
){
    VectorXd pos;
    pos.resize( 3 );
    pos *= 0.0;
    Matrix< double, Eigen::Dynamic, 3 > coord;
    coord.resize( 1, 3 );
    coord( 0, 0 ) = 0.0;
    coord( 0, 1 ) = 0.0;
    coord( 0, 2 ) = 0.0;
    vector< double > charge ( 1, 1.0 );
    return potential_ewald__2D( pos, Tvec, Gvec, coord, charge );
}


//
//
// The ewald-2D method was done with the help of www-user.tu-chemnitz.de/~potts/paper/2dp-ewald.pdf
//
//  "Fast Ewald Summation based on NFFT with Mixed Periodicity" by Franziska Nestler
//
// currently not fast...
//
//

double potential_ewald__2D( 
    VectorXd pos, 
    Matrix3d& Tvec, 
    Matrix3d& Gvec, 
    Matrix< double, Eigen::Dynamic, 3 >& coord, 
    vector< double >& charge
){

   cout << "***************************************************************" << endl;
   cout << " If the 2D ewald method is slow... try changing the alpha value" << endl;
   cout << "***************************************************************" << endl;
 
   if( !suppressoutput ){ 
      cout << "BEGINNING EWALD FOR 2D" << endl;
      cout << "ALPHA : " << alpha << endl;
      cout << fixed << setprecision(10);
   }


   // potEwald   : ewald energy
   // potLRange	: long-ranged part of ewald energy (reciprocal space)
   // potSRange	: short-ranged part of ewald energy (direct space)
   // potSI	: self-interaction part of ewald energy
   // potBK	: background part of ewald energy
   //                - to be used if a uniform background charge needs to
   //                  be subtracted from the input density to obtain the
   //                  ewald energy
   // alpha     : alpha term in ewald summation

   double potEwald = 0.0;
   double potLRange = 0.0;
   double potSRange = 0.0;
   double potSI     = 0.0;
   double potBK     = 0.0;

   //TODO find an alpha that minimized maxG + maxT

   double denom;

   // ... now creating the short-ranged term
   Vector3d T1 = Tvec.row(0);
   Vector3d T2 = Tvec.row(1);
   Vector3d T3 = Tvec.row(2);
   double volume = fabs( T1.dot(T2.cross(T3)) );

      setalpha( Tvec, coord.rows(), volume );

   assert( fabs( T1.dot( T2 ) ) < 1e-13 && "2D EWALD ONLY WORKS FOR PLANAR GEOMETRIES CURRENTLY" );
   assert( fabs( T2.dot( T3 ) ) < 1e-13 && "2D EWALD ONLY WORKS FOR PLANAR GEOMETRIES CURRENTLY" );
   assert( fabs( T1.dot( T3 ) ) < 1e-13 && "2D EWALD ONLY WORKS FOR PLANAR GEOMETRIES CURRENTLY" );
   for( int i = 0; i < coord.rows(); ++i )
      assert( fabs( coord( i, 2 ) ) < 1e-13 && "2D EWALD ONLY WORKS FOR PLANAR GEOMETRIES, SET Z=0" );
   assert( fabs( pos( 2 ) ) < 1e-13 && "2D EWALD ONLY WORKS FOR PLANAR GEOMETRIES, SET Z=0" );
   volume /= sqrt( T3.dot( T3 ) );

   int nPoints = coord.rows();
   double transvecLength;
   Vector3d transvec;
   Vector3d coord1;
   Vector3d coord2;
   Vector3d coordDiff;
   //finding Gvectors to sum over for when argument in exponent is less than 1e-09 = exp(-20)
   Vector3d G1 = Gvec.row(0);
   Vector3d G2 = Gvec.row(1);
   Vector3d G3 = Gvec.row(2);
   Vector3d recipvec, recipvecx, recipvecy;
   double gdotg;
   
   double maxG = 0.0;
   vector<int> mxG (3);
   mxG[ 0 ] = ceil ( alpha * sqrt( volume ) * 4.5 / sqrt( G1.dot( G1 ) ) / M_PI ) + 2;
   mxG[ 1 ] = ceil ( alpha * sqrt( volume ) * 4.5 / sqrt( G2.dot( G2 ) ) / M_PI ) + 2;
//   mxG[ 2 ] = ceil ( sqrt( 20.0 * 4.0 ) * alpha / sqrt( G3.dot( G3 ) ) );
   mxG[ 2 ] = 0; 

   cout << "VOLUME : " << volume << endl;
   cout << "SUMMING OVER : " << mxG[ 0 ] << ", " << mxG[ 1 ] << endl;

   double structureFactor;
   Vector3d relcoord;
   double exparg;
   double rgdotgsqrt;
   for(int ix=-mxG[0];ix<mxG[0]+1;ix++){
      recipvecx = ix * G1;
      for(int iy=-mxG[1];iy<mxG[1]+1;iy++){
         recipvecy = iy * G2;
         for(int iz=-mxG[2];iz<mxG[2]+1;iz++){
            recipvec = recipvecx + recipvecy + iz * G3;
            rgdotgsqrt = sqrt( ix*ix + iy*iy );
            structureFactor *= 0.0;
            if( abs( ix ) + abs( iy ) + abs( iz ) != 0 ){
               for( int iat = 0; iat < nPoints; iat++ ){
                  relcoord(0) = coord( iat, 0 ) - pos(0);
                  relcoord(1) = coord( iat, 1 ) - pos(1);
                  relcoord(2) = coord( iat, 2 ) - pos(2);
                  structureFactor += charge[ iat ] * 1. / 2. / sqrt( volume ) * cos( recipvec.dot( relcoord ) ) * \
                                     ( 2. * erfc( M_PI * rgdotgsqrt / alpha / sqrt( volume ) ) ) \
                                     / rgdotgsqrt; 
                  //cout << structureFactor << " made from " << cos( recipvec.dot( relcoord ) ) << " " << exp( -(gdotg/alpha/alpha/4.0) ) << " " << gdotg << endl;
               }
            }else{
               for( int iat = 0; iat < nPoints; iat++ )
                   structureFactor -= charge[ iat ] * 2. * sqrt( M_PI ) / alpha / sqrt( volume ) / sqrt( volume );
            }
            potLRange += structureFactor ;
         }
      }
   }

   if( !suppressoutput ){
      cout << "LR term " << potLRange << " summed over " << mxG[0] << " " << mxG[1] << " " << mxG[2] << endl;
   }


   // finding tvectors to sum over for when argument in erfc is less than 4.5, giving 1e-10 tolerance 
   // In doing so, we have that the argument inside the erfc() is given by
   //               erfc( alpha * n * T)
   // where n is the number of translations in a given x,y,z direction.  Thus we want to have
   //               alpha * n * T < 4.5
   double maxT = 0.0;
   vector<int> mxT (3);
   mxT[ 0 ] = ceil ( 4.5 / alpha / sqrt( T1.dot( T1 ) ) ) + 2;
   mxT[ 1 ] = ceil ( 4.5 / alpha / sqrt( T2.dot( T2 ) ) ) + 2;
//   mxT[ 2 ] = ceil ( 4.5 / alpha / sqrt( T3.dot( T3 ) ) ); 
   mxT[ 2 ] = 0; 

   for(int ix=-mxT[0];ix<mxT[0]+1;ix++){
      for(int iy=-mxT[1];iy<mxT[1]+1;iy++){
         for(int iz=-mxT[2];iz<mxT[2]+1;iz++){

            for(int idim = 0;idim < 3;idim++){
               transvec(idim) = 1.0 * ix * T1(idim) + 1.0 * iy * T2(idim) + 1.0 * iz * T3(idim); 
            }

            for( int iat = 0; iat < nPoints; iat++ ){
               coord1(0) = coord( iat, 0 ) - pos(0);
               coord1(1) = coord( iat, 1 ) - pos(1);
               coord1(2) = coord( iat, 2 ) - pos(2);
               coordDiff = coord1 + transvec; 
               transvecLength = sqrt( coordDiff.dot( coordDiff ) );
               if( transvecLength > 1e-06 ){
                  potSRange += charge[iat] * (erfc( transvecLength * alpha )/transvecLength);
               }
            }

         }
      }
   }

   //cout << setprecision(8);
   if( !suppressoutput ){
      cout << "SR term " << setw(12) << potSRange << " summed over " << mxT[0] << " " << mxT[1] << " " << mxT[2] << endl;
   }

   //creating neutralizing background charge
   double totcharge = 0.0;
   for ( int iat = 0; iat < nPoints; iat++ ){
      totcharge += charge[ iat ];
   }
   potBK = - 2. * sqrt( M_PI ) / alpha / volume;
   potBK *= 0.0;

   if( !suppressoutput ){
      cout << "BK term " << setw(12) << potBK << endl;
   }
   // seeing whether we are calculating the potential where an ion is located, in which case we need to subtract out
   // a self interaction term from the long ranged potential

   //    first we translate to reference cell
   /* THE FOLLOWING ONLY WORKS FOR ORTHOGONAL LATTICE VECTORS... */
   /*
   for( int i = 0; i < 3; i++ ){
      double tempdotprod = 0.0;
      for( int j = 0; j < 3; j++ ){
         tempdotprod += pos( j ) * Tvec( i, j );
      }
      tempdotprod /= Tvec.row( i ).dot( Tvec.row( i ) );
      if( tempdotprod > 1.0 - 1e-04 ){
         pos( 0 ) -= int( tempdotprod ) * Tvec( i, 0 );
         pos( 1 ) -= int( tempdotprod ) * Tvec( i, 1 );
         pos( 2 ) -= int( tempdotprod ) * Tvec( i, 2 );
      } 
      if( tempdotprod < 0.0 + 1e-04 ){
         pos( 0 ) += int( tempdotprod ) * Tvec( i, 0 );
         pos( 1 ) += int( tempdotprod ) * Tvec( i, 1 );
         pos( 2 ) += int( tempdotprod ) * Tvec( i, 2 );
      } 
   }
   */

   //    see if any atom has that same coordinate 
   for( int iat = 0; iat < nPoints; iat++ ){
      coord1(0) = coord( iat, 0 ) - pos(0);
      coord1(1) = coord( iat, 1 ) - pos(1);
      coord1(2) = coord( iat, 2 ) - pos(2);
      if ( coord1.dot( coord1 ) < 1e-10 ){
         potSI += charge[iat];
      }
   }
   potSI *= - ( 2.0 ) * alpha/sqrt(M_PI);


   if( !suppressoutput ){
      cout << "SE term " << setw(12) << potSI << endl;
   }

   potEwald =  (potSRange + potLRange + potSI) + potBK;
//   potEwald =  (potSRange + potLRange + potSI);
   if( !suppressoutput ){
      cout << potEwald << endl;
   }
   return potEwald;
}
