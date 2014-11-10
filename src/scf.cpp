#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include "scf.h"
#include "scf_diis.h"
#include "ao_ints.h"
#include "Eigen/Dense"
#include "cellinfo.h"

using namespace std;
using namespace Eigen;



void SCF::Initial_Guess( UnitCell& UCell, SuperCell& SCell ){
   //
   //sets up initial density matrix
   //
   bool density_read = false;
   if( pbc.ret_type() == GAMMA ){
      if( readin_density ){
          density_read = read_binary_density( UCell, SCell );
      }
      if( !density_read ) // then we failed to read density or readin_density is not set
        DensXd.irrep(0) = MatrixXd::Identity(nmo_ucell,nmo_ucell);
   }
   else if( pbc.ret_type() == KPOINT ){
      if( readin_density ){
          density_read = read_binary_density( UCell, SCell );
      }
      if( !density_read ) // then we failed to read density or readin_density is not set
        DensXd.irrep(0) = MatrixXd::Identity(nmo_ucell,nmo_ucell);
   }
}


void SCF::SCF_Iteration( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints ){
   clock_t t;
   double energyold=0.0;
   energynew=999.0;
   densdif=999.0;
   evdif = 999.0;

   scf_iter = 0;
   VMatrixXd old_evals = e_vals;
   
   if( readin_mo && read_binary_mo_coeff( UCell, SCell, aoints )){
     nucnuc = aoints.getNucNuc();
     if( pbc.ret_type() == GAMMA  ){ build_dens_gamma( UCell, SCell ); };
     if( pbc.ret_type() == KPOINT ){ 
         build_dens_k( UCell, SCell ); 
         FT_dens_k_to_real( UCell, SCell );
     }
     for( int i = 0; i < nirreps; i++ ){
        OnsiteXd.irrep(i)  = make_onsiteE_matr( UCell, SCell, aoints, i );
        VenXd.irrep(i)  = make_coulomb_matr( UCell, SCell, aoints, ZeroDensXd, i );
        TmatXd.irrep(i) = make_kinetic_matr( UCell, SCell, aoints, i ); 
        CmatXd.irrep(i)  = make_coulomb_matr( UCell, SCell, aoints, DensXd, i );
        EmatXd.irrep(i)  = make_exchange_matr( UCell, SCell, aoints, DensXd, i );
        FockXd.irrep(i) = CmatXd.irrep(i) + EmatXd.irrep(i) + TmatXd.irrep(i) + OnsiteXd.irrep(i);
     }
     return;
   }
   readin_mo = false; // we failed to read in MO coefficients

   cout << "Starting iterations..." << endl;

   while( (scf_iter< max_scf_iter ) && ((fabs(energyold-energynew)>eTol) || (densdif>dTol) || (evdif>evTol)) ){

      energyold = energynew;
      if( scf_iter == 0 ){
         cout << "making nuclear...." << flush;
         nucnuc = aoints.getNucNuc();
         cout << "done." << endl;
         cout << "making onsite, kinetic and ven...." << flush;
         for( int i = 0; i < nirreps; i++ ){
            OnsiteXd.irrep(i)  = make_onsiteE_matr( UCell, SCell, aoints, i );
            VenXd.irrep(i)  = make_coulomb_matr( UCell, SCell, aoints, ZeroDensXd, i );
            TmatXd.irrep(i) = make_kinetic_matr( UCell, SCell, aoints, i ); 
         }
         cout << "done." << endl;
      }

      for(int i=0;i<nirreps;i++){
         CmatXd.irrep(i)  = make_coulomb_matr( UCell, SCell, aoints, DensXd, i );
         EmatXd.irrep(i)  = make_exchange_matr( UCell, SCell, aoints, DensXd, i );
         FockXd.irrep(i) = CmatXd.irrep(i) + EmatXd.irrep(i) + TmatXd.irrep(i) + OnsiteXd.irrep(i);

         old_DensXd.irrep(i) = DensXd.irrep(i);
      }

      // useful for the DIIS
      Create_Gamma_MatrixXd( UCell, SCell, FockXd, SCell_FockXd );

      double damping_factor = 0.05;

      if( pbc.ret_type() == GAMMA ){
         old_evals.irrep(0) = e_vals.irrep(0);
         // need to create SCell_FockXd; by filling it in with the FockXd's
         if ( scf_iter > 0 && use_diis && densdif < 0.10 ){
            MatrixXd fock_with_diis = MatrixXd::Zero( nmo_scell, nmo_scell ) ;
            DIIS.use_diis( &fock_with_diis, SCell_FockXd.irrep(0), SCell_DensXd.irrep(0) );
            SelfAdjointEigenSolver<MatrixXd> eigensolver( fock_with_diis );
            e_vals.irrep(0) = eigensolver.eigenvalues();
            e_vecsXd.irrep(0) = eigensolver.eigenvectors();
         }else{
            if( use_diis ) DIIS.restart();
            SelfAdjointEigenSolver<MatrixXd> eigensolver( SCell_FockXd.irrep(0) );
            e_vals.irrep(0) = eigensolver.eigenvalues();
            e_vecsXd.irrep(0) = eigensolver.eigenvectors();
         }
         build_dens_gamma( UCell, SCell );
         if( densdif >= 0.05 ){
           for( int i = 0; i < nirreps; ++i ){
             DensXd.irrep( i ) = damping_factor * DensXd.irrep( i ) + ( 1. - damping_factor ) * old_DensXd.irrep( i );
           }
         }
      } 


      if( pbc.ret_type() == KPOINT ){

         if ( scf_iter > 0 && use_diis ){
            MatrixXd fock_with_diis = MatrixXd::Zero( nmo_scell, nmo_scell ) ;
            DIIS.use_diis( &fock_with_diis, SCell_FockXd.irrep(0), SCell_DensXd.irrep(0) );
            //copying the gamma point fock matrix to the irrep matrices
            Create_Block_MatrixXd( UCell, SCell, FockXd, fock_with_diis );
            //fourier transforming fock matrices to get the fock matrices in K-space
            FT_fock_real_to_k( UCell, SCell );
            //diagonalizing the k-Focks at each K-point 
            for(int i=0;i<nirreps;i++){
               old_evals.irrep(i) = e_vals.irrep(i);
               SelfAdjointEigenSolver<MatrixXcd> eigensolver(FockXcd.irrep(i));
               e_vals.irrep(i) = eigensolver.eigenvalues();
               e_vecsXcd.irrep(i) = eigensolver.eigenvectors();
            }
         }else{
            if( use_diis ) DIIS.restart();
            //fourier transforming fock matrices to get the fock matrices in K-space
            FT_fock_real_to_k( UCell, SCell );
            //diagonalizing the k-Focks at each K-point 
            for(int i=0;i<nirreps;i++){
               old_evals.irrep(i) = e_vals.irrep(i);
               SelfAdjointEigenSolver<MatrixXcd> eigensolver(FockXcd.irrep(i));
               e_vals.irrep(i) = eigensolver.eigenvalues();
               e_vecsXcd.irrep(i) = eigensolver.eigenvectors();
            }
         }

         //building the density matrix in k-space
         build_dens_k( UCell, SCell );
         
         //doing an inverse fourier transform to get real space density
         FT_dens_k_to_real( UCell, SCell ); 
         Create_Gamma_MatrixXd( UCell, SCell, DensXd, SCell_DensXd );
      } 

      // This forces the density diagonal to be equal to 1 ...
      if( enforce_charge_sym ){
        cout << "Enforcing charge symmetry by setting diagonal elements of Density Matrix to 1..." << endl;
        double charge_max_diff = 0.0;
        for( int i = 0; i < nmo_ucell; ++i ){
          if( fabs( DensXd.irrep( 0 )( i, i ) - 1. ) > charge_max_diff ){
            charge_max_diff = fabs( DensXd.irrep( 0 )( i, i ) - 1. );
          }
          DensXd.irrep( 0 )( i, i ) = 1.;
        }
        printf( "   - MAX DIFFERENCE : %24.16e \n", charge_max_diff );
      }
      // Enforcing inversion symmetry to get rid of pesky wrapping errors ...
      Enforce_Inversion_Symmetry( UCell, SCell, aoints, DensXd );
      Create_Gamma_MatrixXd( UCell, SCell, DensXd, SCell_DensXd );

      densdif = 0.0;
      for(int i = 0; i < nirreps; i++){
         densdif+=((old_DensXd.irrep(i) - DensXd.irrep(i)).lpNorm<Infinity>());
      }

      evdif=0.0;
      evdif+=((old_evals.irrep(0) - e_vals.irrep(0)).lpNorm<Infinity>());



      energynew=0.0;
      for(int i=0;i<nirreps;i++){
         for(int j=0;j<nmo_ucell;j++){
            for(int k=0;k<nmo_ucell;k++){
               energynew += (DensXd.irrep(i)(j,k) * (FockXd.irrep(i)(j,k) + TmatXd.irrep(i)(j,k) + VenXd.irrep(i)(j,k) + OnsiteXd.irrep(i)(j,k)));
            }
         }
      }
      energynew *= 0.5;
      energynew /= nirreps;

      energydif = fabs(energyold - energynew);
      printf("SCF ITERATION: %3d  ENERGY: %20.14e  dE: %20.14e  dDens: %20.14e  dEig: %20.14e \n",scf_iter,energynew,energydif,densdif, evdif);
      fflush( stdout );
      E_Energy = energynew;
     
      if( write_matrices ){
        printPartialMatrix( "VEN", VenXd );
        printPartialMatrix( "KINETIC", TmatXd );
        printPartialMatrix( "ONSITE", OnsiteXd );
        printPartialMatrix( "COULOMB", CmatXd );
        printPartialMatrix( "EXCHANGE", EmatXd );
        printPartialMatrix( "FOCK", FockXd );
        printPartialMatrix( "NEW DENSITY", DensXd );
      }
      scf_iter++;
   }

   if( scf_iter == max_scf_iter ){ printf( "SCF FAILED : CONVERGENCE NOT REACHED IN %3d CYCLES! EXITING... \n", max_scf_iter ); exit( EXIT_FAILURE ); };

   // constructing fock from the current density
   for(int i=0;i<nirreps;i++){
      CmatXd.irrep(i)  = make_coulomb_matr( UCell, SCell, aoints, DensXd, i );
      EmatXd.irrep(i)  = make_exchange_matr( UCell, SCell, aoints, DensXd, i );
      FockXd.irrep(i) = CmatXd.irrep(i) + EmatXd.irrep(i) + TmatXd.irrep(i) + OnsiteXd.irrep(i);
   }

   write_binary_density( UCell, SCell );
   write_binary_mo_coeff( UCell, SCell, aoints );

}



void SCF::get_Total_Energy(){
   one_body_energy = 0.0;
   two_body_energy = 0.0;
   total_dipole_correction = 0.0;
   K_Energy = 0.0;
   EN_Energy = 0.0;
   EE_Energy = 0.0;
   XC_Energy = 0.0;
   Onsite_Energy = 0.0;
   energynew=0.0;
   /* shitty implementation, but easy to see if correct */
   for(int i=0;i<nirreps;i++){
      for(int j=0;j<nmo_ucell;j++){
         for(int k=0;k<nmo_ucell;k++){
            Onsite_Energy += (DensXd.irrep(i)(j,k)) * ( OnsiteXd.irrep(i)(j,k) );
            K_Energy += (DensXd.irrep(i)(j,k)) * ( TmatXd.irrep(i)(j,k) );
            EN_Energy += (DensXd.irrep(i)(j,k)) * ( VenXd.irrep(i)(j,k) );
            XC_Energy += (DensXd.irrep(i)(j,k)) * ( EmatXd.irrep(i)(j,k) );
            EE_Energy += (DensXd.irrep(i)(j,k)) * ( CmatXd.irrep(i)(j,k) - VenXd.irrep( i )( j, k ) );
            energynew += (DensXd.irrep(i)(j,k) * (FockXd.irrep(i)(j,k) + TmatXd.irrep(i)(j,k) + VenXd.irrep(i)(j,k) + OnsiteXd.irrep(i)(j,k)));
         }
      }
   }
   energynew *= 0.5;
   XC_Energy *= 0.5;
   EE_Energy = 0.5 * EE_Energy;
   one_body_energy = EN_Energy + K_Energy + Onsite_Energy;
   two_body_energy = EE_Energy + XC_Energy;

   ground_state_energy = energynew + nucnuc;
   printf("::::::::: SCF END :::::::::::\n" );
   printf("FINAL ENERGY                : %20.16f \n", ground_state_energy );
   printf("FINAL ENERGY + DIPOLE CORR. : %20.16f \n", ( ground_state_energy + total_dipole_correction ) );
   printf("ONSITE ENERGY               : %20.16f \n", Onsite_Energy );
   printf("HARTREE ENERGY              : %20.16f \n", EE_Energy );
   printf("EXCHANGE ENERGY             : %20.16f \n", XC_Energy );
   printf("NUCLEAR ATT. ENERGY         : %20.16f \n", EN_Energy );
   printf("KINETIC ENERGY              : %20.16f \n", K_Energy  );
   printf("NUCLEAR-NUCLEAR ENERGY      : %20.16f \n", nucnuc );
   printf("DIPOLE CORRECTION           : %20.16f \n", total_dipole_correction );
   printf("TOTAL ELECTROSTATIC CONT.   : %20.16f \n", ( EN_Energy + EE_Energy + nucnuc ) );
   printf("TOTAL ELECT. MINUS DIPOLE   : %20.16f \n", ( EN_Energy + EE_Energy + nucnuc - total_dipole_correction ) );
   printf("-------------------------------------\n");
   printf("TOTAL ONE-BODY CONTRIBUTION : %20.16f \n", one_body_energy );
   printf("TOTAL TWO-BODY CONTRIBUTION : %20.16f \n", two_body_energy );
   printf("TOTAL ZER-BODY CONTRIBUTION : %20.16f \n", nucnuc );
   printf("SUM                         : %20.16f \n", (one_body_energy + nucnuc + two_body_energy ) );
   if( pbc.ret_type() == GAMMA ){
      printf("EIGENVALUES           :        \n");
      for( int i = 0; i < nmo_scell; ++i ){
         printf("%5d %20.16f \n", i, e_vals.irrep( 0 )(i, 0 ) );
      }
      cout << setprecision( 16 ) << fixed;
      cout << "BAND GAP (au,eV) : " << ( e_vals.irrep( 0 )( nocc, 0 ) - e_vals.irrep( 0 )( nocc - 1, 0 ) ) << " " << \
                                       ( e_vals.irrep( 0 )( nocc, 0 ) - e_vals.irrep( 0 )( nocc - 1, 0 ) ) * 27.211 << endl;
   }
   if( pbc.ret_type() == KPOINT ){
      for( int j = 0; j < nirreps; ++j ){
        printf("EIGENVALUES AT BAND %d :        \n", j );
        for( int i = 0; i < nmo_ucell; ++i )
            printf("%5d %14.8f \n", i, e_vals.irrep( j )(i, 0 ) );
      } 
      cout << endl;
   }
}


void SCF::Compute_Energy( UnitCell &UCell, SuperCell& SCell, aoIntegralFactory& aoints ){

   /* Initializes matrices, gets SCF options, and sets up SCF to use a specific model (a la PPP) 
    * if necessary, also gets integrals possibly, and sets up for K-Point or Gamma calculation
    */
   Init( UCell, SCell );

   // Creates initial density matrix 
   Initial_Guess( UCell, SCell );
   
   // Goes through iterations of the SCF process until stop criteria are met
   SCF_Iteration( UCell, SCell, aoints ); 

   //Possibly check MO orthogonality?
  
   //Possibly check if k-space and r-space energy are equal for periodic systems

   // Gets total energy if the various convergence criteria were met (including checking MO orthogonality
   // possibly, calling to the nuclear energy that is dependant on the model
   get_Total_Energy();

   //Checks energy using MO integrals
//   check_energy();

   //Write out the MO integrals
//   write_out_MOs_ewald( "MOintegrals.dat" );
}

