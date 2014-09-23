#ifndef OPTIONS_H
#define OPTIONS_H

#include <map>
#include <vector>
#include <string>

struct Options {
    void    AddOptionDouble( const char* option, double      default_val );
    void    AddOptionInt   ( const char* option, int         default_val );
    void    AddOptionBool  ( const char* option, bool        default_val );
    void    AddOptionString( const char* option, const char* default_val );
    void    AddOptionStringWF( const char* option, void (*func)( const char* ), const char* default_val );

    void    SetOptionDouble( const char* option, double      val );
    void    SetOptionInt   ( const char* option, int         val );
    void    SetOptionBool  ( const char* option, const char* val );
    void    SetOptionBool  ( const char* option, bool        val );
    void    SetOptionString( const char* option, const char* val );

    void    SetOptionsFromList( std::vector< std::string > str );

    void    PrintOptions();

    double      GetOptionDouble( const char* option );
    int         GetOptionInt   ( const char* option );
    bool        GetOptionBool  ( const char* option );
    const char* GetOptionString( const char* option );

    std::map< const char*, double      > option_double_map;
    std::map< const char*, int         > option_int_map;
    std::map< const char*, bool        > option_bool_map;
    std::map< const char*, const char* > option_string_map;

    std::map< const char*, double      >::iterator iter_double;
    std::map< const char*, int         >::iterator iter_int;
    std::map< const char*, bool        >::iterator iter_bool;
    std::map< const char*, const char* >::iterator iter_string;

};

#endif
