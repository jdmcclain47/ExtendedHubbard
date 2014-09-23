#ifndef SCF_H
#define SCF_H

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include "cellinfo.h"
#include "common.h"
#include <cstring>
#include "scf_diis.h"
#include "ao_ints.h"
#include "options.h"
#include "Eigen/Dense"
// include number of irreps/noocc and w/e with some other options

//typedef Eigen::Matrix<Eigen::MatrixXd,Eigen::Dynamic,1> VMatrixXd;

class SCF {
    public:
      void Compute_Energy( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints );

      bool read_binary_density( UnitCell& UCell, SuperCell& SCell, const char* density_file );
      bool read_binary_density( UnitCell& UCell, SuperCell& SCell );
      void write_binary_density( UnitCell& UCell, SuperCell& SCell, const char* density_file );
      void write_binary_density( UnitCell& UCell, SuperCell& SCell );

      bool  read_binary_mo_coeff( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory&, const char* in_file );
      bool  read_binary_mo_coeff( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory&  );
      void write_binary_mo_coeff( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory&, const char* in_file );
      void write_binary_mo_coeff( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& );

      int get_nocc(){return nocc;};
      int get_nmo_scell(){return nmo_scell;};
      double get_rcut(){ return rcut; };
      double get_ground_state_energy(){ return ground_state_energy; };
      double get_kinetic_energy(){ return K_Energy; };
      double get_xc_energy(){ return XC_Energy; };
      double get_onsite_energy(){ return Onsite_Energy; };
      double get_one_body_energy(){ return one_body_energy; };
      double get_two_body_energy(){ return two_body_energy; };
      double get_nuclear_energy(){ return nucnuc; };
      bool get_converge_pseudo_1d(){ return converge_pseudo_1d; };
      bool get_converge_pseudo_2d(){ return converge_pseudo_2d; };
      double get_ewald_se(){ return ewald_self_interaction; };
      double get_eTol(){ return eTol; };
      VMatrixXd& get_gamma_evecs(){return e_vecsXd;};
      VMatrixXcd& get_kpoint_evecs(){return e_vecsXcd;};
      VMatrixXd& get_evals(){return e_vals;};
      PBC_CLASS& get_pbc(){ return pbc; };
   private:
      void printPartialMatrix( const char* title, VMatrixXd& in_matr );
      void Init( UnitCell& UCell, SuperCell& SCell );
      void Initial_Guess( UnitCell& UCell, SuperCell& SCell );
      void SCF_Iteration( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints );
      void get_Total_Energy(); 

      void setDefaultOptions();
      Options Opts;

      int ndiis;
      diis DIIS;
      PBC_CLASS pbc; 
      double E_Energy;
      double nucnuc;
      double N_Energy;
      double T_Energy;
      double K_Energy, XC_Energy, EE_Energy, EN_Energy, Onsite_Energy;
      double ground_state_energy, one_body_energy, two_body_energy, total_dipole_correction;
      double densdif, evdif;
      double dTol, eTol, evTol;
      int max_scf_iter;
      double energydif;
      double energynew;
      int scf_iter;
      int nocc;
      int nmo_ucell;   //number of morbitals in unit cell
      int nmo_scell;   //number of morbitals in super cell
      int nirreps;

      double itol, rcut;
      bool ws_cut, ws_cut_correction, use_ws;
      bool converge_pseudo_1d;
      bool converge_pseudo_2d;
      bool sph_cut, sph_cut_correction;
      bool readin_density, readin_mo, read_density_diagonal;
      bool do_gamma, do_kpoint;
      bool use_diis;
      bool do_xc;
      bool write_matrices;
      bool xcewald;
      bool add_coulomb_correction;
      double ewald_self_interaction;

      /* A list of matrices, the number of which is the number of irreducibles
       */
      VMatrixXd CmatXd;   //Coulomb Matrix
      VMatrixXd SmatXd;   //Overlap Matrix
      VMatrixXd EmatXd;   //Exchange Matrix
      VMatrixXd TmatXd;   //Kinetic Matrix
      VMatrixXd OnsiteXd;   //Onsite Matrix for PPP
      VMatrixXd FockXd;   //Fock Matrix
      VMatrixXd VenXd;   //nuclear attraction matrix, added to find the electron energy
      VMatrixXd DensXd;   //Density Matrix
      VMatrixXd ZeroDensXd;   //Density Matrix
      VMatrixXd old_DensXd;   //Density Matrix at previous iteration

      VMatrixXd e_vals;
      VMatrixXd e_vecsXd;
      VMatrixXcd e_vecsXcd;

/*
      Eigen::MatrixXd atomic_dist_matr;
      Eigen::MatrixXi ws_matr;
      Eigen::Matrix3d inv_matr;
      Eigen::MatrixXd ewald_pot_matr;
      Eigen::MatrixXi scell_matr;
      Eigen::MatrixXd kernel_matr;
*/

      Eigen::MatrixXd make_exchange_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, VMatrixXd& inDens, int which_unit_cell );
      Eigen::MatrixXd make_coulomb_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, VMatrixXd& inDens, int which_unit_cell );
      Eigen::MatrixXd make_kinetic_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, int which_unit_cell );
      Eigen::MatrixXd make_onsiteE_matr( UnitCell& UCell, SuperCell& SCell, aoIntegralFactory& aoints, int which_unit_cell );

      /* The corresponding complex versions of the fock and density
       * along with the fock matrix for the entire supercell
       */
      VMatrixXcd FockXcd;
      VMatrixXcd DensXcd;

      /* The Gamma point matrices encompassing the entire supercell
       * in one matrix
       */
      VMatrixXd SCell_FockXd;
      VMatrixXd SCell_DensXd;


      void check_energy();
      void FT_fock_real_to_k( UnitCell& UCell, SuperCell& SCell );
      void Create_Gamma_MatrixXd( UnitCell& UCell, SuperCell& SCell, VMatrixXd& block_matrix, VMatrixXd& gamma_matrix );
      void Create_Block_MatrixXd( UnitCell& UCell, SuperCell& SCell, VMatrixXd& block_matrix, VMatrixXd& gamma_matrix );
      void Create_Block_MatrixXd( UnitCell& UCell, SuperCell& SCell, VMatrixXd& block_matrix, Eigen::MatrixXd& gamma_matrix );
      void build_dens_k( UnitCell& UCell, SuperCell& SCell );
      void build_dens_gamma( UnitCell& UCell, SuperCell& SCell );
      void FT_dens_k_to_real( UnitCell& UCell, SuperCell& SCell );

      void write_out_density();
      void read_density();
      void write_out_MOs_ewald( const char* filename );

};

#endif
