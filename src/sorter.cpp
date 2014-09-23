#include <vector>
#include <utility>
#include <algorithm>

static bool pair_compare( std::pair< size_t, double > pair1, std::pair< size_t, double > pair2 ){
    return ( pair1.first < pair2.first );
}

void sorter(
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    const size_t& arr_size
){
    std::pair< size_t, double > indx_pair;
    std::vector< std::pair< size_t, double > > indx_pair_arr;
    indx_pair_arr.resize( arr_size );
    for( size_t i = 0; i < arr_size; ++i ){
      indx_pair = std::make_pair( indx_arr[ i ], dble_arr[ i ] );
      indx_pair_arr[ i ] = indx_pair;
    } 
    std::sort( indx_pair_arr.begin(), indx_pair_arr.begin() + arr_size, pair_compare );
    for( size_t i = 0; i < arr_size; ++i ){
      indx_arr[ i ] = indx_pair_arr[ i ].first;
      dble_arr[ i ] = indx_pair_arr[ i ].second;
    } 
}
