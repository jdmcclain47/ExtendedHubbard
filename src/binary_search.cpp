#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

int binary_search( 
    std::vector< size_t >& indx_arr,
    size_t INDEX_TO_BE_FOUND,
    int lower_bound,
    int higher_bound
){
  if( lower_bound > higher_bound ) return -1;
  size_t mid = ( lower_bound + higher_bound ) / 2;
  if( INDEX_TO_BE_FOUND < indx_arr[ mid ] )
    return binary_search( indx_arr, INDEX_TO_BE_FOUND, lower_bound, mid-1 );
  else if( INDEX_TO_BE_FOUND > indx_arr[ mid ] )
    return binary_search( indx_arr, INDEX_TO_BE_FOUND, mid+1, higher_bound );
  else
    return mid;
}
