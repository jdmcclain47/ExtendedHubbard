#ifndef CREATE_HEADER_H
#define CREATE_HEADER_H

#include <iostream>
#include <stdio.h>
#include "ao_ints.h"
#include "common.h"
#include "cellinfo.h"
#include <cstring>

std::string create_kernel_suffix(
  UnitCell& UCell,
  SuperCell& SCell,
  aoIntegralFactory& aoints
);

char* make_kernel_header(
  UnitCell& UCell,
  SuperCell& SCell,
  aoIntegralFactory& aoints,
  size_t& size
);

char* read_kernel_header();

#endif
