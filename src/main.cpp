#include <iostream>
#include "binary_search.h"
#include <iomanip>
#include "parser.h"
#include "ao_ints.h"
#include "davidson_new.h"
#include "mp2.h"
#include "cellinfo.h"
#include "create_header.h"
#include "math.h"
#include "scf.h"
#include "davidson.h"
#include <stdio.h>
#include <vector>
#include <string>
#include "common.h"
#include "mointegrals.h"
#include "moint_read_file.h"
#include <time.h>
#include <complex>
#include "moint_kpoint_twobody.h"

#include "mpi.h"

int main( int argc, char* argv[] ){

   int numtasks, taskid, source;
   MPI_Init( &argc, &argv );
   MPI_Comm_size( MPI_COMM_WORLD, &numtasks );
   MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
   if( taskid == 0 )
     std::cout << "RUNNING " << numtasks << " TASKS WITH MPI." << std::endl;
   std::cout << "TASK " << taskid << " STARTED!" << std::endl;

   UnitCell            UCell;
   SuperCell           SCell;
   aoIntegralFactory   aoints;
   SCF                 scf;
   Davidson            cis;
   moIntegralFactory moints;

   

   if( taskid == 0 ){
    
     UCell.Init();
     SCell.Init( UCell );
     aoints.Init( UCell, SCell, "PARAM", "PARAM" );
    
     scf.Compute_Energy( UCell, SCell, aoints );

     //aoints.set_non_coulomb_kernel_to_zero();

   }

   if( ( scf.get_pbc() ).type == KPOINT ){
     if( taskid == 0 ){
       cis.Init( UCell, SCell, moints );
       moint_driver_k kmoints;
       kmoints.kpoint_init( UCell, SCell, aoints, scf.get_evals(), scf.get_kpoint_evecs() );
       std::vector< std::complex< double > > kintegrals = kmoints.create_kpoint_2e_ints_fast_FULL( UCell, SCell, aoints, scf.get_kpoint_evecs() );
       kmoints.check_exchange_energy( UCell, SCell, aoints, kintegrals );
       printf( "  - EXCHANGE ENERGY FROM SCF      : %20.16f \n", scf.get_xc_energy() );
       kmoints.write_out_moints( kintegrals, ".KMOINTS" );
       cis.set_evals_and_moints( scf.get_evals(), kintegrals );
       cis.DavidsonCIS_kpoint( 30000, -6 );
     }
   } else {
     clock_t t;
     if( taskid == 0 ){
       moints.Init( UCell, SCell, aoints, scf.get_pbc(), scf.get_evals(), scf.get_gamma_evecs(), scf.get_kpoint_evecs() );
       t=clock();
       moints.write_gamma_moint_to_fcidump( "FCIDUMP", -10 );
       t=clock()-t; std::cout << "MOINTS NEW TIME " << (double)t/CLOCKS_PER_SEC << std::endl;
     }
     std::string mointfile ( ".MOINTEGRALS" );



     double mp2energy;
     cis.Init( UCell, SCell, moints );

     //
     //
     // for full gamma mo integrals
     //
     //

     //if( taskid == 0 ){
     //  t=clock();
     //  moints.write_gamma_moint_to_file( ".MOINTS" );
     //  t=clock()-t; std::cout << "MOINTS NEW TIME " << (double)t/CLOCKS_PER_SEC << std::endl;

     //  std::vector< double > allints = moints.read_gamma_moints_from_file( ".MOINTS" );
     //  //moints.check_gamma_exchange_energy( allints, scf.get_xc_energy() );
     //  //mp2energy = mp2_gamma_full_moints( UCell, SCell, scf.get_evals(), allints );
     //  //printf( "MP2 ENERGY          : %20.16f \n", mp2energy );
     //  //printf( "MP2 ENERGY D KPOINT : %20.16f \n", mp2energy / SCell.nkpt );
     //  cis.set_evals_and_moints( scf.get_evals(), allints );
     //  cis.DavidsonCIS( 3000, 4, 10, mointfile );
     //}
     //cis.CheckerCIS();
//Davidson_Gamma_MPI( UCell, SCell, aoints, scf.get_gamma_evecs(), scf.get_evals(), 4, 10 );
MP2_Gamma_MPI( UCell, SCell, aoints, scf.get_gamma_evecs(), scf.get_evals(), 10 );
MPI_Finalize();
return;

     //
     //
     // for index p gamma mo integrals 
     //
     //

     //t=clock();
     //moints.write_gamma_mointb_ind_p( mointfile, -10 );
     //t=clock()-t; std::cout << "MO NEW NEW TIME " << (double)t/CLOCKS_PER_SEC << std::endl;
     //clock_t mp2timer;
     //if( taskid == 0 ) mp2timer = clock();
     //mp2energy = mp2_gamma_ind_p( UCell, SCell, scf.get_evals(), mointfile, false, moints, argc, argv );
     //if( taskid == 0 ){
     //  mp2timer = clock() - mp2timer;
     //  printf( "MP2 ENERGY    : %20.16f \n", mp2energy );
     //  printf( "MP2 EN / NKPT : %20.16f \n", mp2energy / SCell.nkpt );
     //  printf( "MP2 TIME      : %20.16f \n", mp2timer / (double)CLOCKS_PER_SEC );
     //}
     //cis.set_evals( scf.get_evals() );
     //cis.DavidsonCIS( 3000, 1, 10, mointfile );
     //mp2energy = mp2_gamma_ind_p( UCell, SCell, scf.get_evals(), mointfile, true, moints );
     //printf( "MP2 ENERGY  : %20.16f \n", mp2energy );
     //printf( "MP2 ENERGY D KPOINT : %20.16f \n", mp2energy / SCell.nkpt );


     //
     //
     // making FCIDUMP file for use with ACES3
     //
     //
     if( taskid == 0 ){
       t=clock();
       moints.write_gamma_moint_to_fcidump_binary( "FCIDUMP_BINARY", -10 );
       t=clock()-t; std::cout << "MOINTS NEW TIME " << (double)t/CLOCKS_PER_SEC << std::endl;
     }
   }

   MPI_Finalize();

}
