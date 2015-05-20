#ifndef AOINTS_CLASS_H
#define AOINTS_CLASS_H

#include <vector>
#include "cellinfo.h"
#include "Eigen/Dense"
#include "options.h"
#include "ppp.h"
#include <cstring>

enum pppkernel { 
  MN, // Mataga-Nishimoto
  COULOMB, // 1/r
  OHNO, // Ohno-Klopman
  MAZUMDAR // Ohno-Klopman like with a dielectric constant
};

enum kernel { 
  WS, // Wigner-Seitz
  EWALD, // ewald kernel with no G=0 correction
  EWALDMSE, // ewald kernel minus self energy
  NONE // set equal to zero
};

struct inversionSym {
    std::vector< int > invij;
    int naoUnitCell, naoSuperCell; 

    void makeUniqueList( std::vector< Eigen::Vector3i >& tList, int nx, int ny, int nz );
    void setup( UnitCell& UCell, SuperCell& SCell );
};

class aoIntegralFactory {
    inversionSym inversion;
    int dim;
    int naoUnitCell, naoSuperCell; 
    int naoUnit_sqp1;
    int nTrans;
    int tolEwald;
    bool do_ppp_kernel_correction, use_new_ints;
    bool do_xc_ppp_correction, do_coulomb_ppp_correction;
    bool printMatr, readChkpt, printMatrToFile;
    double pppkerneldist;
    double ewald_self;
    double nucnuc;
    double prefac;
    Eigen::Matrix3d invmat;
    std::vector< double > aoDistMatr;    
    std::vector< double > aoHopMatr;
    std::vector< double > aoEwaldMatr;
    std::vector< double > aoNonCoulombMatr;
    std::vector< double > aoCutMatrWS;
    Eigen::MatrixXi scell_matr;
    pppkernel pppkern;
    kernel xckern; 
    kernel coulombkern; 
    kernel corr1kern, corr2kern; 
    Options Opts;
    PPPModel PPP;

    void setDefaultOptions();
    void makeHoppingMatr( UnitCell& UCell, SuperCell& SCell );
    void make_scell_matr( SuperCell& SCell );
    void makeDistMatr( UnitCell& UCell, SuperCell& SCell, bool explicit_MIC );
    void makeEwaldMatr( UnitCell& UCell, SuperCell& SCell, int tol );
    void makeCutMatrWS( UnitCell& UCell, SuperCell& SCell );
    void setKernel( kernel& inkern, const char* kernName );
    int getPerElement( int iao, int jao, int which_cell );

    /* PPP begin (more in public) */

    void SetPPPKernel( std::string inkern ){ SetPPPKernel( inkern.c_str() ); };
    void SetPPPKernel( pppkernel inkern ){ pppkern = inkern; };
    pppkernel GetPPPKernel()             { return pppkern; };

    std::string getKernelStr( const kernel& inkern );

    void SetXCKernel( std::string inkern ){ SetXCKernel( inkern.c_str() ); };
    void SetXCKernel( kernel inkern ){ xckern = inkern; };
    kernel GetXCKernel()             { return xckern; };

    void SetCoulombKernel( std::string inkern ){ SetCoulombKernel( inkern.c_str() ); };
    void SetCoulombKernel( kernel inkern ){ coulombkern = inkern; };
    kernel GetCoulombKernel()             { return coulombkern; };

    double getPPPInt( std::string ao1, std::string ao2, double dist );
    double getPPPInt( std::string ao1, std::string ao2, double dist, const pppkernel inkern );
    void makeNonCoulombMatr(
      UnitCell& UCell,
      SuperCell& SCell,
      double rcut
    );

    bool read_ewald( UnitCell& UCell, SuperCell& SCell, int tolEwald, std::string filename );
    void write_ewald( UnitCell& UCell, SuperCell& SCell, int tolEwald, std::string filename );

    /* PPP end */

    public :

      void Init( UnitCell& UCell, SuperCell& SCell, const char* ao_options_filename, const char* ppp_options_filename );
      void printPerSuperMatrix( SuperCell& SCell, std::vector< double >& inmatr, const char* title );
      void printPerSuperMatrix( SuperCell& SCell, std::vector< double >& inmatr, const char* title, int width, int prec );
      void printCorr1ToFile( SuperCell& SCell, const char* title );
      void printCorr2ToFile( SuperCell& SCell, const char* title );
      bool get_integral_proc(){ return use_new_ints; };
      int getInvTrans( int which_cell );
      void set_non_coulomb_kernel_to_zero(){ for( int i = 0; i < aoNonCoulombMatr.size(); ++i ){ aoNonCoulombMatr[ i ] = 0.0;} };
 
      double getCoulombInt( int iat, int jat, int which_cell );
      double getXCInt( int iat, int jat, int which_cell );
      double getCorr1Int( int iat, int jat, int which_cell );
      double getCorr2Int( int iat, int jat, int which_cell );
      double getSelfInt(){ return ewald_self; };
      double getHopping( int iat, int jat, int which_cell ){ return aoHopMatr[ getPerElement( iat, jat, which_cell ) ]; };
      double getNucNuc(){ return nucnuc; };

      /* PPP begin */

      std::string getPPPKernelStr();
      std::string getXCKernelStr(){ return getKernelStr( xckern ); };
      std::string getCoulombKernelStr(){ return getKernelStr( coulombkern ); };
      void SetPPPKernel( const char* inkern );
      void SetXCKernel( const char* inkern );
      void SetCorr1Kernel( const char* inkern );
      void SetCorr2Kernel( const char* inkern );
      void SetCoulombKernel( const char* inkern );
      double getHubbardU( std::string ao1 ){ return PPP.get_hubbard_u( ao1 ); };
      double getOnsiteE( std::string ao1 ){ return PPP.get_hubbard_e( ao1 ); };

      /* PPP end */

};

#endif
