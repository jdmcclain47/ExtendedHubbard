#include "options.h"
#include <string>
#include <string.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <map>


void Options::AddOptionDouble( const char* option, double default_val ){
    option_double_map[ option ] = default_val;
}

void Options::AddOptionInt( const char* option, int default_val ){
    option_int_map[ option ] = default_val;
}

void Options::AddOptionBool( const char* option, bool default_val ){
    option_bool_map[ option ] = default_val;
    //const char* tf = ( default_val ) ? "TRUE" : "FALSE";
    //printf( "Adding bool option '%s' - default set to '%s' \n", option, tf );
}

void Options::AddOptionString( const char* option, const char* default_val ){
    option_string_map[ option ] = default_val;
}

void Options::AddOptionStringWF( const char* option, void (*func)( const char* ), const char* default_val ){
    option_string_map[ option ] = default_val;
    (*func)( default_val );
}

void Options::SetOptionDouble( const char* option, double      val ){
    iter_double = option_double_map.find( option );    
    if( iter_double == option_double_map.end() ){ // we didn't find the option
      printf( "WARNING : Option '%s' not found for SetOptionDouble! ", option );
      printf( "Should first use AddOptionDouble( '%s', default_value ). \n", option );
    }
    option_double_map[ option ] = val;
}

void Options::SetOptionInt   ( const char* option, int         val ){
    iter_int = option_int_map.find( option );    
    if( iter_int == option_int_map.end() ){ // we didn't find the option
      printf( "WARNING : Option '%s' not found for SetOptionDouble! ", option );
      printf( "Should first use AddOptionDouble( '%s', default_value ). \n", option );
    }
    option_int_map[ option ] = val;
}

void Options::SetOptionBool   ( const char* option, const char*      val ){
    iter_bool = option_bool_map.find( option );    
    if( iter_bool == option_bool_map.end() ){ // we didn't find the option
      printf( "WARNING : Option '%s' not found for SetOptionBool! ", option );
      printf( "Should first use AddOptionBool( '%s', default_value ). \n", option );
    }
    if( strcmp( val, "TRUE" ) == 0 ){
        option_bool_map[ option ] = true; return;
    }
    if( strcmp( val, "FALSE" ) == 0 ){
        option_bool_map[ option ] = false; return;
    }
    printf( "ERROR : For SetOptionBool, option '%s' expects either 'TRUE' or 'FALSE'! \n", option );
    exit( EXIT_FAILURE );
}

void Options::SetOptionBool   ( const char* option, bool      val ){
    iter_bool = option_bool_map.find( option );    
    if( iter_bool == option_bool_map.end() ){ // we didn't find the option
      printf( "WARNING : Option '%s' not found for SetOptionBool! ", option );
      printf( "Should first use AddOptionBool( '%s', default_value ). \n", option );
    }
    option_bool_map[ option ] = val; 
}

void Options::SetOptionString( const char* option, const char* val ){
    iter_string = option_string_map.find( option );    
    if( iter_string == option_string_map.end() ){ // we didn't find the option
      printf( "WARNING : Option '%s' not found for SetOptionDouble! ", option );
      printf( "Should first use AddOptionDouble( '%s', default_value ). \n", option );
    }
    option_string_map[ option ] = val;
}

double Options::GetOptionDouble( const char* option ){
        iter_double =  option_double_map.find( option );    
    if( iter_double == option_double_map.end() ){ // we didn't find the option
      printf( "ERROR : Option '%s' not found! \n", option );
      exit( EXIT_FAILURE );
    }
    return iter_double->second;
}

int Options::GetOptionInt   ( const char* option ){
        iter_int =  option_int_map.find( option );    
    if( iter_int == option_int_map.end() ){ // we didn't find the option
      printf( "ERROR : Option '%s' not found! \n", option );
      exit( EXIT_FAILURE );
    }
    return iter_int->second;
}

bool Options::GetOptionBool( const char* option ){
        iter_bool =  option_bool_map.find( option );    
    if( iter_bool == option_bool_map.end() ){ // we didn't find the option
      printf( "ERROR : Option '%s' not found! \n", option );
      exit( EXIT_FAILURE );
    }
    return iter_bool->second;
}

const char* Options::GetOptionString( const char* option ){
        iter_string =  option_string_map.find( option );    
    if( iter_string == option_string_map.end() ){ // we didn't find the option
      printf( "ERROR : Option '%s' not found! \n", option );
      exit( EXIT_FAILURE );
    }
    return iter_string->second;
}

void Options::PrintOptions(){
    printf( "-------------------------------------------------\n" );
    printf( "OPTIONS\n" );
    printf( "-------------------------------------------------\n" );
    for( iter_double =  option_double_map.begin(); 
         iter_double != option_double_map.end(); 
       ++iter_double )
    {
      printf( "%25s %10c %20.16f \n", iter_double->first, ' ', iter_double->second );
    }
    for( iter_int =  option_int_map.begin(); 
         iter_int != option_int_map.end(); 
       ++iter_int )
    {
      printf( "%25s %10c %20d \n", iter_int->first, ' ', iter_int->second );
    }
    for( iter_bool =  option_bool_map.begin(); 
         iter_bool != option_bool_map.end(); 
       ++iter_bool )
    {
      const char* tf = ( iter_bool->second ) ? "TRUE" : "FALSE";
      printf( "%25s %10c %20s \n", iter_bool->first, ' ', tf );
    }
    for( iter_string =  option_string_map.begin(); 
         iter_string != option_string_map.end(); 
       ++iter_string )
    {
      printf( "%25s %10c %20s \n", iter_string->first, ' ', iter_string->second );
    }
    printf( "-------------------------------------------------\n" );
    printf( "\n" );
}

void Options::SetOptionsFromList( std::vector< std::string > str ){

    // searching this list for any options taking type double
    for( iter_double =  option_double_map.begin(); 
         iter_double != option_double_map.end(); 
       ++iter_double )
    {
      for( int i = 0; i < str.size(); ++i ){
        if( strcmp( iter_double->first, str[ i ].c_str() ) == 0 ){ // we found the option!
          double num;
          std::istringstream( str[ i + 1 ] ) >> num;
          SetOptionDouble( iter_double->first, num );
        }
      }
    }

    // searching this list for any options taking type int 
    for( iter_int =  option_int_map.begin(); 
         iter_int != option_int_map.end(); 
       ++iter_int )
    {
      for( int i = 0; i < str.size(); ++i ){
        if( strcmp( iter_int->first, str[ i ].c_str() ) == 0 ){ // we found the option!
          int num;
          std::istringstream( str[ i + 1 ] ) >> num;
          SetOptionInt( iter_int->first, num );
          break;
        }
      }
    }

    // searching this list for any options taking type bool 
    for( iter_bool =  option_bool_map.begin(); 
         iter_bool != option_bool_map.end(); 
       ++iter_bool )
    {
      for( int i = 0; i < str.size(); ++i ){
        if( strcmp( iter_bool->first, str[ i ].c_str() ) == 0 ){ // we found the option!
          SetOptionBool( iter_bool->first, str[ i + 1 ].c_str() );
          break;
        }
      }
    }

    // searching this list for any options taking type string 
    for( iter_string =  option_string_map.begin(); 
         iter_string != option_string_map.end(); 
       ++iter_string )
    {
      for( int i = 0; i < str.size(); ++i ){
        if( strcmp( iter_string->first, str[ i ].c_str() ) == 0 ){ // we found the option!
          SetOptionString( iter_string->first, str[ i + 1 ].c_str() );
          break;
        }
      }
    }
}
