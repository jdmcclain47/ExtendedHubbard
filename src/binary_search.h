#ifndef BINARY_SEARCH_H
#define BINARY_SEARCH_H

#include <vector>
#include <math.h>

int binary_search( 
    std::vector< size_t >& indx_arr,
    size_t INDEX_TO_BE_FOUND,
    int lower_bound,
    int higher_bound
);

#endif
