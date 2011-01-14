//---------------------------------------------------------------------------
//
// itl wrappers, callable from C, C++, and Fortran
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// Copyright Notice
// + 2010 University of Chicago
//
//--------------------------------------------------------------------------
#include "itl.h"

// C and C++ wrapper
void ITL_begin() {

}

// fortran wrapper
extern "C"
void itl_begin_() {

  fprintf(stderr, "Hello world from ITL\n");

}
//--------------------------------------------------------------------------
