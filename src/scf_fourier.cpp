#include <vector>
#include <complex>
#include <iostream>
#include "scf.h"
#include "common.h"
#include "Eigen/Dense"
#include "cellinfo.h"
using namespace std;
using namespace Eigen;

void SCF::FT_fock_real_to_k( UnitCell& UCell, SuperCell& SCell ){
   MatrixXcd temp_matr;
   temp_matr.resize(nmo_ucell,nmo_ucell);
   for(int i=0;i<nirreps;i++){
      /* the following is to create the fock in k-space at the k-point labelled by kpoint(i)
       */
      temp_matr *= 0.0;
      for(int j=0;j<nmo_ucell;j++){
         for(int k=0;k<nmo_ucell;k++){
            for(int rn=0;rn<SCell.total_number_of_cells;rn++){
               real(temp_matr(j,k)) += cos( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * FockXd.irrep(rn)(j,k);
               imag(temp_matr(j,k)) += sin( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * FockXd.irrep(rn)(j,k);
            }
         }
      }
      FockXcd.irrep(i) = temp_matr;
   }
}


void SCF::FT_dens_k_to_real( UnitCell& UCell, SuperCell& SCell ){
   MatrixXcd temp_matr;
   temp_matr.resize(nmo_ucell,nmo_ucell);
   for(int rn=0;rn<SCell.total_number_of_cells;rn++){
      /* the following is to create the density matrix in r-space from the one in k-space 
       */
      temp_matr *= 0.0;
      for(int j=0;j<nmo_ucell;j++){
         for(int k=0;k<nmo_ucell;k++){
            for(int i=0;i<nirreps;i++){
               real(temp_matr(j,k)) += cos( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * real(DensXcd.irrep(i)(j,k));
               real(temp_matr(j,k)) += sin( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * imag(DensXcd.irrep(i)(j,k));

               imag(temp_matr(j,k)) += cos( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * imag(DensXcd.irrep(i)(j,k));
               imag(temp_matr(j,k)) -= sin( (SCell.kpoints[i]).dot(SCell.translations[rn]) ) * real(DensXcd.irrep(i)(j,k));
            }
            if(fabs(imag(temp_matr(j,k))) > 1e-10){
               cout << "***********  WARNING  **********************" << endl;
               cout << "Real density matrix has imaginary component." << endl;
               cout << "At Rn = "<< rn << "(" << j << "," << k << "): " << temp_matr(j,k) << endl;
               cout << " .. This may be due to SOMO states .. " << endl << endl;
            }
            DensXd.irrep(rn)(j,k) = real(temp_matr(j,k)) / ( 1. * nirreps );
         }
      }
   }
}
