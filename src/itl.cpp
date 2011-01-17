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

  int rank; // MPI usual
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    fprintf(stderr, "Beginning ITL\n");

}

// C and C++ wrapper
void ITL_get_data() {

}

// fortran wrapper
extern "C"
void itl_get_data_(int *nelv, int *nx1, int *ny1, int *nz1, double *vx, double *vy, double *vz, double *xm1, double *ym1, double *zm1) {

  int rank; // MPI usual
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {

    fprintf(stderr, "nelv = %d nx1 = %d ny1 = %d nz1 = %d\n", *nelv, *nx1, *ny1, *nz1);

    fprintf(stderr, "velocity field:\n");
    for (int i = 0; i < *nx1 * *ny1 * *nz1; i++) {
	  fprintf(stderr, "[vx vy vz] = %.3lf %.3lf %.3lf at [x y z] = %.3lf %.3lf %.3lf\n", 
		  vx[i], vy[i], vz[i], 
		  xm1[i], ym1[i], zm1[i]);
    }

  }

}
//--------------------------------------------------------------------------
