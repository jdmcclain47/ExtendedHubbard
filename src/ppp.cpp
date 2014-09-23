#include <vector>
#include "parser.h"
#include "ppp.h"
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <cstring>
#include <cassert>
#include <iostream>
#include <map>
#include <math.h>
#include <algorithm>
#include <sstream>

static double au2ev = 27.211396132;

void PPPModel::Init( const char* infilename ){
    filename = infilename;
    std::vector< std::string > params = ParseText( filename );
    Init( params, filename );
}

void PPPModel::Init( std::vector< std::string >& params, const char* infilename ){
    filename = infilename;
    double dtemp;
    std::vector< std::string >::iterator params_iter;
    std::vector< std::string >::iterator temp_iter;
    std::string req_string ( "NAO" );
    params_iter = find( params.begin(), params.end(), req_string );
    if( params_iter == params.end() ){
      printf( "ERROR : For the function PPPModel::Init, the filename '%s' "
              "does not have the required string '%s' \n",
              filename, req_string.c_str() );
    }
    params_iter++; // moving on to the next iterator that should be the number of AO's
    std::istringstream( *params_iter ) >> nao;
    params_iter++; // moving on to the next iterator that should be the list of AO's
    for( int i = 0; i < nao; ++i ){
      aolist.push_back( *params_iter ); 
      params_iter++;
    } 

    /* onsite U elements */
    std::string uString ( "PPP_U" );
    onsiteU.resize( nao );
    setValues1VD( uString, &onsiteU, 0.0, params );
    for( int i = 0; i < nao; ++i ) onsiteU[ i ] /= au2ev; // converting from eV to au

    /* onsite E elements */ 
    std::string eString ( "PPP_E" );
    onsiteE.resize( nao );
    setValues1VD( eString, &onsiteE, 0.0, params );
    for( int i = 0; i < nao; ++i ) onsiteE[ i ] /= au2ev; // converting from eV to au

    /* prefac for coulomb interactions */ 
    std::string prefacString ( "PREFAC" );
    setValueDouble( prefacString, &prefac, 1.0, params );

    /* tolerance for hopping integral */ 
    std::string ttolString ( "TTOL" );
    setValueDouble( ttolString, &ttol, -1.0, params );
    dtol = pow( 10., ttol );

    /* hopping matrix elements */
    std::string thop ( "THOP" );
    setValuesHopping( thop, &hopping, params );
    for( hop_iter = hopping.begin(); hop_iter != hopping.end(); ++hop_iter ){ // converting from eV to au
      (hop_iter->second).val /= au2ev;
    }

    /* dielectric constant */
    std::string dielString ( "DIEL" );
    setValueDouble( dielString, &diel, 2.0, params );

    print();
}

int PPPModel::getOrbLocation( std::string aoName ){
    ao_iter = find( aolist.begin(), aolist.end(), aoName );
    return (ao_iter - aolist.begin() );
}

int PPPModel::getOrbLocation( std::string aoName1, std::string aoName2 ){
    ao_iter = find( aolist.begin(), aolist.end(), aoName1 );
    int n1 = ( ao_iter - aolist.begin() );
    ao_iter = find( aolist.begin(), aolist.end(), aoName2 );
    int n2 = ( ao_iter - aolist.begin() );
    if( n1 > n2 )
      return ( (n1*(n1+1))/2+n2 );
    else
      return ( (n2*(n2+1))/2+n1 );
}

double PPPModel::getValuesOnsiteU( std::string aoName ){
    int iLoc = getOrbLocation( aoName );
    return onsiteU[ iLoc ];
}

double PPPModel::getValuesOnsiteE( std::string aoName ){
    int iLoc = getOrbLocation( aoName );
    return onsiteE[ iLoc ];
}

double PPPModel::getValuesHopping( std::string ao1, std::string ao2, double distance ){
    int iLoc = getOrbLocation( ao1, ao2 );
    for( hop_iter = hopping.equal_range( iLoc ).first; hop_iter != hopping.equal_range( iLoc ).second; ++hop_iter )
      if( fabs( ( (*hop_iter).second ).dist - distance ) < dtol ) return ( ( (*hop_iter).second ).val );
    return 0.0;
}

void PPPModel::getOrbString( int iKey, std::string* orb1, std::string* orb2 ){
    int ival = 0;
    int ivalsq = 0;
    int itot = 0;
    bool found = false;
    while( !found ){
      ival++;
      if( itot + ival > iKey ){
        found = true;
        (*orb2) = aolist[ iKey - itot ];
        (*orb1) = aolist[ ival - 1 ];
      }else{
        itot += ival;
      } 
    }
}

void PPPModel::setValueDouble( std::string varName, double* val, double default_val, std::vector< std::string >& params ){
    std::vector< std::string >::iterator params_iter;
    params_iter = find( params.begin(), params.end(), varName );
    if( params_iter == params.end() ){
      printf( "WARNING : For the function PPPModel::setValueDouble, no variable "
              "'%s' found in filename '%s'.\n"
              "          Default set to : %14.f \n", 
              varName.c_str(), filename, default_val );    
      (*val) = default_val;
    }
    else{
      params_iter++;
      std::istringstream( *params_iter ) >> (*val);
    }
}

void PPPModel::setValues1VD( std::string varName, std::vector< double >* varContainer, double default_val, std::vector< std::string >& params ){
    int nVal = (*varContainer).size();
    int iLoc;
    double dtemp;
    bool found = false;
    for( int i = 0; i < nVal; ++i ) (*varContainer)[ i ] = default_val;
    std::vector< std::string >::iterator params_iter;
    std::vector< std::string >::iterator temp_iter;

    std::vector< std::string >::iterator params_current = params.begin();
    while( params_current != params.end() ){
      params_iter = find( params_current, params.end(), varName );
      if( params_iter != params.end() ){
        params_iter++; // looking at next element that should be the ao name
        std::string aoLookup = *params_iter;
        temp_iter = find( aolist.begin(), aolist.end(), aoLookup );
        iLoc = (temp_iter-aolist.begin() );
        if( temp_iter == aolist.end() ){
          printf( "ERROR : For the function PPPModel::Init (filename '%s') "
                  "in setting '%s', no orbital '%s' found. \n", 
                  filename, varName.c_str(), aoLookup.c_str() );
          exit( EXIT_FAILURE );
        }
        params_iter++; // looking at next element that should be the value 
        if( params_iter == params.end() ){
          printf( "ERROR : For the function PPPModel::Init (filename '%s') "
                  "in setting '%s' for orbital '%s' EOF reached. \n",
                  filename, varName.c_str(), aoLookup.c_str() );
          exit( EXIT_FAILURE );
        }
        std::istringstream( *params_iter ) >> dtemp;
        (*varContainer)[ iLoc ] = dtemp;
        found = true;
      }
      params_current = params_iter;
    }

    if( !found ){
      printf( "WARNING : For the function PPPModel::setValues1VD, no variable "
              "'%s' found in filename '%s'.\n"
              "          Default set to : %14.f \n",
              varName.c_str(), filename, default_val );    
    }
}


void PPPModel::setValuesHopping( 
    std::string varName, 
    std::multimap< int, hopping_pair >* varContainer, 
    std::vector< std::string >& params 
){
    int iLoc;
    double dtemp;
    bool found = false;
    std::vector< std::string >::iterator params_iter;
    std::vector< std::string >::iterator temp_iter;

    std::vector< std::string >::iterator params_current = params.begin();
    while( params_current != params.end() ){
      params_iter = find( params_current, params.end(), varName ); // finding the variable name
      if( params_iter != params.end() ){
        params_iter++; // looking at next element that should be the ao name
        std::string aoLookup1 = *params_iter;
        params_iter++; // looking at next element that should be the ao name
        std::string aoLookup2 = *params_iter;
        iLoc = getOrbLocation( aoLookup1, aoLookup2 ); // assigning a key to this orbital pair

        hopping_pair thop;
        params_iter++; // looking at next element that should be the distance; 
        assert( params_iter != params.end() && "ERROR : Hopping elements should be in "
                "the form ATOM1 ATOM2 DISTANCE VALUE" );
        std::istringstream( *params_iter ) >> thop.dist;
        params_iter++; // looking at next element that should be the hopping value; 
        assert( params_iter != params.end() && "ERROR : Hopping elements should be in "
               "the form ATOM1 ATOM2 DISTANCE VALUE" );
        std::istringstream( *params_iter ) >> thop.val;

        (*varContainer).insert( std::make_pair( iLoc, thop ) ); // inserting into container
        found = true;
      }
      params_current = params_iter;
    }

    if( !found ){
      printf( "ERROR : No hopping integrals found in filename '%s'. Aborted. \n", filename );
    }

}


void PPPModel::print(){

    std::string atomtitle ( "ORBITAL" );
    std::string valuetitle ( "VALUE (a.u.)" );
    std::string disttitle ( "DISTANCE (a.u.)" );

    printf( "--------------------------------------------\n" );
    printf( "PPP PARAMETERS                              \n" );
    printf( "--------------------------------------------\n" );

    printf( "ONSITE U ELEMENTS \n" );
    printf( "\t %8s %20s \n", atomtitle.c_str(), valuetitle.c_str() ); 
    for( int i = 0; i < nao; ++i )
      printf( "\t %8s %20.16f \n", aolist[ i ].c_str(), onsiteU[ i ] );

    printf( "ONSITE E ELEMENTS \n" );
    printf( "\t %8s %20s \n", atomtitle.c_str(), valuetitle.c_str() ); 
    for( int i = 0; i < nao; ++i )
      printf( "\t %8s %20.16f \n", aolist[ i ].c_str(), onsiteE[ i ] );

    printf( "PREFAC            : %20.16f \n", prefac );
    printf( "HOPPING TOLERANCE : 10^(%f) \n", ttol );

    printf( "HOPPING ELEMENTS \n" );
    printf( "\t %8s %8s %20s %20s \n", atomtitle.c_str(), atomtitle.c_str(), disttitle.c_str(), valuetitle.c_str() ); 
    for( hop_iter = hopping.begin(); hop_iter != hopping.end(); ++hop_iter ){ 
      std::string orb1, orb2;
      getOrbString( hop_iter->first, &orb1, &orb2 ); 
      printf( "\t %8s %8s %20.16f %20.16f \n", orb1.c_str(), orb2.c_str(), (hop_iter->second).dist, (hop_iter->second).val );
    }
    printf( "--------------------------------------------\n\n" );
}
