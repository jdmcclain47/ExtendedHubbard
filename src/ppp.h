#ifndef PPP_H
#define PPP_H

#include <map>
#include <vector>
//#include "parser.h"
//#include <cstdlib>
//#include <stdio.h>
//#include <string>

struct hopping_pair {
    double dist;
    double val;
};

struct PPPModel {

    int nao; // number of unique types of orbitals
    std::vector< std::string > aolist;
    std::vector< std::string >::iterator ao_iter;
    const char* filename;
    double prefac;
    double diel;
    double ttol;
    double dtol;
    std::vector< double > onsiteU; // the U param in the PPP/Hubbard model
    std::vector< double > onsiteE; // a constant energy shift for a particular site
    std::multimap< int, hopping_pair > hopping;
    std::multimap< int, hopping_pair >::iterator hop_iter;

    void Init( const char* infilename );
    void Init( std::vector< std::string >& params, const char* infilename );
    void setValues1VD( std::string varName, std::vector< double >* varContainer, double default_val, std::vector< std::string >& params );
    void setValueDouble( std::string varName, double* val, double default_val, std::vector< std::string >& params );
    int getOrbLocation( std::string aoName );
    int getOrbLocation( std::string aoName1, std::string aoName2 );
    void getOrbString( int iKey, std::string* orb1, std::string* orb2 );
    void setValuesHopping( 
        std::string varName, 
        std::multimap< int, hopping_pair >* varContainer, 
        std::vector< std::string >& params 
    );
    double getValuesHopping( std::string ao1, std::string ao2, double distance );
    double getValuesOnsiteU( std::string aoName );
    double getValuesOnsiteE( std::string aoName );

    double get_diel(){ return diel; };
    double get_hopping( std::string ao1, std::string ao2, double distance ){ return getValuesHopping( ao1, ao2, distance ); };
    double get_hubbard_u( std::string aoName ){ return getValuesOnsiteU( aoName ); };
    double get_hubbard_e( std::string aoName ){ return getValuesOnsiteE( aoName ); };
    double get_prefac(){ return prefac; };
    void print();

};

#endif
