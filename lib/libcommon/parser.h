#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <string>

std::vector< std::string > ParseText( const char * infile );
std::vector< std::string > ParseText( const char * infile, bool quiet_mode );

#endif
