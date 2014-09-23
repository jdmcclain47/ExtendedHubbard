#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstring>
#include "parser.h"
using namespace std;

vector< std::string > ParseText( const char* infile ){
    return ParseText( infile, 0 );
}

vector< std::string > ParseText( const char * infile, bool quiet_mode )
{

    vector< string > strvec;
    string line, temp, str;
    int i, j, whichline;
    size_t maxln;
    ifstream inread;
    inread.open( infile );
    if( inread.is_open() ){
        cout << "Parsing text file '" << infile << "'..." << endl;
    }else{
        cout << "ERROR : Error opening text file '" << infile << "'..." << endl;
        exit( EXIT_FAILURE );
    }

    whichline = 0;
    while( getline( inread, line ) )
    {
    if( line.length() > 0 ) // ... ignore our line if it's empty
    {
        // ... removing all comments
        maxln = line.find( "#" );
        cout << "line #" << whichline << " (len="<< line.length() << "): " << line << endl;
        if( maxln == -1 ) // ... we have no comment!
        {
           i = 0;
           j = 0;
           while( j < line.length() )
           {
             while( isspace( line[ i ] ) ) // ... eliminating trailing zeros before a string
             {
                 i++;
             } 
             if( i == line.length() ) break; // ... end our search if the beginning of the next string is the
                                             //     end of the line
             j = i;
             while( !isspace( line[ j ] ) && j < line.length() ) // ... while the location of the end of the string is
                                                                 //     less than the line length and it is not an empty
                                                                 //     space, we're still looking at a string
             {
                 j++; 
             }
             str = line.substr( i, (j-i) );
             if( str.length() == 1 && str[ 0 ] == '-' ){
               cout << "ERROR : Error parsing file '" << infile << "' at line " << whichline << ". " <<
                       "Stray '-' found.\n" << 
                       "Either remove the '-' or remove any spaces between it and the next string." << endl; 
               exit( EXIT_FAILURE );
             }
             strvec.push_back( str ); 
             //cout << line.substr( i, (j-i) ) << endl; // ... string starts at position 'i' and has length 'j-i'
             i = j; // ... if we aren't at the end of the line, we begin our search again, starting 
                    //     at the whitespace
           }
            
        }
        else // ... we do have a comment!
        {
           // ... taking the part of the string before this comment
           cout << "   - ignoring comment : '" << line.substr( maxln, line.length() - maxln ) << "'." << endl;
           line = line.substr( 0, maxln ); 
           // ... trimming all leading spaces
           i = 0;
           j = 0;
           while( j < line.length() )
           {
             while( isspace( line[ i ] ) ) // ... eliminating trailing zeros before a string
             {
                 i++;
             } 
             if( i == line.length() ) break; // ... end our search if the beginning of the next string is the
                                             //     end of the line
             j = i;
             while( !isspace( line[ j ] ) && j < line.length() ) // ... while the location of the end of the string is
                                                                 //     less than the line length and it is not an empty
                                                                 //     space, we're still looking at a string
             {
                 j++; 
             }
             str = line.substr( i, (j-i) );
             if( str.length() == 1 && str[ 0 ] == '-' ){
               cout << "ERROR : Error parsing file '" << infile << "' at line " << whichline << ". " <<
                       "Stray '-' found.\n" << 
                       "Either remove the '-' or remove any spaces between it and the next string." << endl; 
               exit( EXIT_FAILURE );
             }
             strvec.push_back( str );
             //cout << line.substr( i, (j-i) ) << endl; // ... string starts at position 'i' and has length 'j-i'
             i = j; // ... if we aren't at the end of the line, we begin our search again, starting 
                    //     at the whitespace
           }
        }
    }
    else
    {
        cout << "line #" << whichline << " (EMPTY): " << endl;
    }
    whichline++;
    }
    if( inread.eof() ){
        cout << "File read successfully!" << endl;
    }else{
        cout << "ERROR : Error reading file '" << infile << "'. Didn't read to EOF." << endl;
        exit( EXIT_FAILURE );
    }
    inread.close();

    return strvec;
}

