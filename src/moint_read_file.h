#ifndef MOINT_READ_FILE
#define MOINT_READ_FILE

#include <string>
#include <vector>

void read_gamma_mointb_ind_p( 
    std::string ofilename,
    const int& p_index,
    std::vector< double >& dble_arr,
    std::vector< size_t >& indx_arr,
    size_t& size_arr,
    bool& optimize_size_arr
);

#endif
