#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include "scf_diis.h"
#include "Eigen/Dense" 

using namespace std;
using namespace Eigen;

static int diisiter;

void diis::Init( int nmax_, int nblocks_ ){
   nmax = nmax_;
   nblocks = nblocks_;

   diisiter = 1;
   focklist.resize( nmax );
   errlist.resize( nmax );
}

void diis::restart(){
   diisiter = 1;
}

void diis::use_diis( MatrixXd* outFock, MatrixXd inFock, MatrixXd inDens){
   // where we update the diis fock list
   int updatehere = diisiter - 1;
   // number of diis matrices
   int ndiis = diisiter;
   // size of diismatrix (ndiis + 1)
   int diissize;

   if ( diisiter > nmax ){
      updatehere = (diisiter - 1) % nmax;
      ndiis = nmax;
   }
   diissize = ndiis + 1;

   // updating fock list
   focklist[ updatehere ] = inFock;

/*
   for(int i = 0; i < ndiis;i ++ ){
      cout << focklist[i] << endl << endl;
   }
*/

   // now creating the error matrix and updating error list
   int fockdim = inFock.rows();
   MatrixXd errMatrix = MatrixXd::Zero( fockdim, fockdim );
   errMatrix = inDens * inFock - inFock * inDens;
   errlist[ updatehere ] = errMatrix;

   // checking if we have at least 2 fock matrices to perform diis
   if ( ndiis > 1 ){
      // now creating the diis matrix
      MatrixXd diisMatrix = MatrixXd::Zero( diissize, diissize );

      for ( int i = 0; i < diissize; i ++ ){
         for ( int j = 0; j <= i; j ++ ){
            if( i == ( diissize - 1) || j == ( diissize - 1)){
               diisMatrix(i,j) = -1.0;
               diisMatrix(j,i) = -1.0;
            } 
            else{
               double temp = ( errlist[i] * errlist[j] ).trace();
               diisMatrix(i,j) = temp;
               diisMatrix(j,i) = temp;
            }
         }
      }
      diisMatrix( diissize - 1, diissize - 1) = 0.0;

      // now solving for coefficients
      VectorXd rhsvector = VectorXd::Zero( diissize );
      rhsvector( diissize - 1 ) = -1;
      VectorXd coeffs;
      coeffs = diisMatrix.colPivHouseholderQr().solve( rhsvector );

      // now returning fock matrix
      for ( int i = 0; i < fockdim; i ++ ){
         for ( int j = 0; j <= i ; j ++ ){
            (*outFock) (i,j) = 0.0;
            (*outFock) (j,i) = 0.0;
         }  
      }
      // do over ndiis, not diissize
      for ( int i = 0; i < ndiis; i++){
         (*outFock) += coeffs( i ) * focklist[ i ];
      }
   }
   else{
      (*outFock) = inFock; 
   }
   diisiter ++;

}
/*
int main(){
   diis DIIS;
   DIIS.Init(2, 1);

   MatrixXd f1 = MatrixXd::Random(5,5);
   MatrixXd f2 = MatrixXd::Random(5,5);
   MatrixXd f3 = MatrixXd::Random(5,5);
   MatrixXd f4 = MatrixXd::Random(5,5);
   MatrixXd d1 = MatrixXd::Random(5,5);
   MatrixXd d2 = MatrixXd::Random(5,5);
   MatrixXd d3 = MatrixXd::Random(5,5);
   MatrixXd d4 = MatrixXd::Random(5,5);
   
   MatrixXd asdf = MatrixXd::Zero(5,5);
   DIIS.use_diis( &asdf, f1, d1);
   cout << asdf << endl;
   DIIS.use_diis( &asdf, f2, d2);
   cout << asdf << endl;
   DIIS.use_diis( &asdf, f3, d3);
   cout << asdf << endl;
   DIIS.use_diis( &asdf, f4, d4);
   cout << asdf << endl;
}
*/
