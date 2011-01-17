/*---------------------------------------------------------------------------
 *
 * itl C, C++ interface
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * Copyright Notice
 * + 2010 University of Chicago
 *
--------------------------------------------------------------------------*/

#ifndef _ITL
#define _ITL

#include "mpi.h"

/*-------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void ITL_begin();
#ifdef __cplusplus
extern "C"
#endif
void ITL_get_data();

/*-------------------------------------------------------------------------*/

#endif
