/**
 * @file ItlNekMoab.cpp
 * ITL wrappers, callable from C, C++, and Fortran
 * Created on: June 20, 2011
 * @author Abon
 * @author Tom Peterka
 */

//#define BYTE_SWAP

#include <mpi.h>

#include "diy.h"
#include "util.hpp"
#include "assignment.hpp"
#include "blocking.hpp"
#include "io.hpp"
#include "merge.hpp"

#include "bil.h"

#include "ITL_header.h"
#include "ITL_base.h"

#include <vector>

#include "iMesh.h"
#include "MBiMesh.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "moab/ParallelComm.hpp"

// vtk headers
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkDoubleArray.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkUnstructuredGrid.h"

using namespace std;
using namespace moab;

// vortex
struct vortex_t {
  double x, y, z; 	// location
  int t; 		// time
  double l2; 		// lambda-2 value
};

static int nb; 			// number of local blocks
static int nv; 			// number of local vertices in all local blocks
static int rank = 0; 		// MPI rank
static int groupsize = 1; 	// MPI communicator size

static double *x, *y, *z; 	// geometry
static double **sc; 		// scalar fields
static double **vx, **vy, **vz; // vector fields, with separate components
static int nsf, nvf; 		// number of scalar and vector fields

static vector<vortex_t> vortices;

//--------------------------------------------------------------------------
//
// initialize CIAN
//
// nblocks: number of local blocks in this process
// num_verts_per_block: number of vertices per block
// num_scalar_fields: number of scalar fields for which to allocate space
//  indexed by scalar_id, 0-num_scalar_fields - 1
// num_vector_fields: number of vector fields for which to allocate space

// C and C++ wrapper
void ITL_begin() {

}

// fortran wrapper
extern "C"
void itl_begin_(int *num_blocks, int *num_verts_per_block,
		 int *num_scalar_fields, int *num_vector_fields) {

  nb = *num_blocks;
  nv = nb * *num_verts_per_block;
  nsf = *num_scalar_fields;
  nvf = *num_vector_fields;
  x = new double[nv];
  y = new double[nv];
  z = new double[nv];
  sc = new double*[nsf];
  vx = new double*[nvf];
  vy = new double*[nvf];
  vz = new double*[nvf];
  for (int i = 0; i < nsf; i++)
    sc[i] = new double[nv];
  for (int i = 0; i < nvf; i++) {
    vx[i] = new double[nv];
    vy[i] = new double[nv];
    vz[i] = new double[nv];
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &groupsize);

  if (rank == 0)
    fprintf(stderr, "Beginning CIAN. Number of blocks = %d "
	    "Number of vertices per block = %d Number of scalar fields = %d "
	    "Number of vector fields = %d\n", 
	    nb, *num_verts_per_block, nsf, nvf);

}
//--------------------------------------------------------------------------
//
// terminate ITL
//

// C and C++ wrapper
void ITL_end() {

}

// fortran wrapper
extern "C"
void itl_end_() {

  if (rank == 0)
    fprintf(stderr, "Ending CIAN\n");

}

//--------------------------------------------------------------------------
// 
// get scalar data
//
// scalar_id: index of scalar field [0 - num_scalar_fields - 1]
//  eg. temperature = index 0, pressure = index 1, etc.
// vals: array of scalar values
// nvals: number of values
// ofst: starting offset (index) of values
//

// C and C++ wrapper
void ITL_scalar_data() {

}

// fortran wrapper
extern "C"
void itl_scalar_data_(int *scalar_id, double *vals, int *nvals, int *ofst) {

  for (int i = 0; i < *nvals; i++)
    sc[*scalar_id][*ofst + i] = vals[*ofst + i];

  // only printing rank 0 to reduce the amount of output
//   if (rank == 0) {

//     fprintf(stderr, "CIAN: scalar field:\n");
//     for (int i = 0; i < *nvals; i++) {
// //     for (int i = 0; i < 100; i++) // only printing first 100 vals
//       if (fabs(vals[*ofst + i]) > 2)
// 	  fprintf(stderr, "val[%d] = %.3e\n", i, vals[*ofst + i]);
//     }

//   }

}
//--------------------------------------------------------------------------
// 
// get vector data
//
// vector_id: index of vector field [0 - num_vector_fields - 1]
//  eg. velocity = index 0, flux = index 1, etc.
// x_vals, y_vals, z_vals: arrays of vector components
// nvals: number of values in a single arrays of one component
// ofst: starting offset (index) of values
// derived_scalar_id: index of derived scalar field
//  -1 : don't derive anything
//  >= 0 : derive a scalar, magnitude for now, and save in this scalar field
//  ensure that enough scalar fields were defined in CIAN_begin
//

// C and C++ wrapper
void ITL_vector_data() {

}

// fortran wrapper
extern "C"
void itl_vector_data_(int *vector_id, double *x_vals, double *y_vals, 
		       double *z_vals, int *nvals, int *ofst, 
		       int *derived_scalar_id) {

  for (int i = 0; i < *nvals; i++) {

    vx[*vector_id][*ofst + i] = x_vals[*ofst + i];
    vy[*vector_id][*ofst + i] = y_vals[*ofst + i];
    vz[*vector_id][*ofst + i] = z_vals[*ofst + i];

    // derived field, magnitude for now
    if (*derived_scalar_id >= 0) {
      sc[*derived_scalar_id][*ofst + i] =
	sqrt(vx[*vector_id][*ofst + i] * vx[*vector_id][*ofst + i] +
	     vy[*vector_id][*ofst + i] * vy[*vector_id][*ofst + i] +
	     vy[*vector_id][*ofst + i] * vz[*vector_id][*ofst + i]);
    }
  }

  // debug
  if (rank == 0) {

    fprintf(stderr, "CIAN: vector field:\n");
    for (int i = 0; i < 64; i++)
      fprintf(stderr, "vx, vy, vz = %.3e, %.3e, %.3e\n", 
	      vx[*vector_id][*ofst + i], vy[*vector_id][*ofst + i], 
	      vz[*vector_id][*ofst + i]);

  }

}

//--------------------------------------------------------------------------
// 
// get moab imesh handle
//
// mesh: the moab mesh instance

// C and C++ wrapper
void ITL_mesh() {

}

// fortran wrapper
extern "C"
void itl_mesh_(iMesh_Instance *mesh) {

  int ierr;

  fprintf(stderr, "CIAN_mesh\n");


  // only printing rank 0 to reduce the amount of output
  if (rank == 0) {

    for (int dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
      int numd;
      iMesh_getNumOfType(*mesh, 0, dim, &numd, &ierr);
      fprintf(stderr, "Number of %dd elements = %d\n", dim, numd);
    }

  }

  // cast the imesh instance to a moab instance
  moab::Interface *mb = ((MBiMesh *)*mesh)->mbImpl;

  // get the block sizes (dims)
  Tag tag_x, tag_y, tag_z, tag_dims;
  EntityHandle root_set = mb->get_root_set();
  int dims[3];
  ErrorCode merr;

  merr = mb->tag_get_handle("SEM_DIMS", 3, MB_TYPE_INTEGER, tag_dims);
  assert(merr == MB_SUCCESS);

  merr = mb->tag_get_data(tag_dims, &root_set, 1, dims);
  assert(merr == MB_SUCCESS);

  // get the tag handles for the x, y, z geometry arrays
  int tot_pts = dims[0] * dims[1] * dims[2];

  merr = mb->tag_get_handle("SEM_X", tot_pts, MB_TYPE_DOUBLE, tag_x);
  assert(merr == MB_SUCCESS);

  merr = mb->tag_get_handle("SEM_Y", tot_pts, MB_TYPE_DOUBLE, tag_y);
  assert(merr == MB_SUCCESS);

  merr = mb->tag_get_handle("SEM_Z", tot_pts, MB_TYPE_DOUBLE, tag_z);
  assert(merr == MB_SUCCESS);

  // get blocks
  Range blocks;
  merr = mb->get_entities_by_type(0, MBHEX, blocks);
  assert(merr == MB_SUCCESS);

  // get x,y,z vertex geometry for each block
  int count;
  double *x, *y, *z;

  mb->tag_iterate(tag_x, blocks.begin(), blocks.end(), count, (void *&)x);
  assert(merr == MB_SUCCESS);
  assert(count == (int)blocks.size());

  mb->tag_iterate(tag_y, blocks.begin(), blocks.end(), count, (void *&)y);
  assert(merr == MB_SUCCESS);
  assert(count == (int)blocks.size());

  mb->tag_iterate(tag_z, blocks.begin(), blocks.end(), count, (void *&)z);
  assert(merr == MB_SUCCESS);
  assert(count == (int)blocks.size());

  // print vertex geometry
  fprintf(stderr, "CIAN: MOAB geometry %d blocks\n", count);
//   for (int i = 0; i < count; i++) {
  for (int i = 0; i < 1; i++) { // print only 1 block for now
    for (int j = 0; j < tot_pts; j++)
      fprintf(stderr, "%.3lf %.3lf %.3lf\n", x[i * tot_pts + j], y[i * tot_pts + j], 
	      z[i * tot_pts + j]);
  }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//
// renders a scalar field using vtk
//
// scalar_id: index of scalar field [0 - num_scalar_fields - 1]
// nx, ny, nz: block size
//
// C and C++ wrapper
void ITL_vis_scalar() {

}

// fortran wrapper
extern "C"
void itl_vis_scalar_(int *scalar_id, int *nx, int *ny, int *nz) {

  // init vtk float array
  vtkDoubleArray *coords = vtkDoubleArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(nv);
  for (int i = 0; i < nv; i++)
    coords->SetTuple3(i, x[i], y[i], z[i]);

  // init vtk points
  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);

  // init vtk cells
  vtkCellArray *cells = vtkCellArray::New();
  for (int i = 0; i < nb; i++) { // for all local blocks
    int s = i * *nx * *ny * *nz;
    for (int j = 0; j < *nx * *ny * *nz; j++) {
      // skip the last entry in each row, the last row in each face, and the
      // last face of each block
      if ((j % *nx == *nx - 1) || ((j / *nx) % *ny == *ny - 1) ||
	  (*nz > 1 && (j / (*nx * *ny)) % *nz == *nz - 1))
	continue;
      if (*nz > 1) // 3d
	cells->InsertNextCell(8);
      else // 2d
	cells->InsertNextCell(4);
      cells->InsertCellPoint(s + j);
      cells->InsertCellPoint(s + j + 1);
      cells->InsertCellPoint(s + j + *nx + 1);
      cells->InsertCellPoint(s + j + *nx);
      if (*nz > 1) {
	cells->InsertCellPoint(s + j + *nx * *ny);
	cells->InsertCellPoint(s + j + *nx * *ny + 1);
	cells->InsertCellPoint(s + j + *nx + *nx * *ny + 1);
	cells->InsertCellPoint(s + j + *nx + *nx * *ny);
      }
    }

  } // local blocks

  // debug
  fprintf(stderr, "vtk scalar field\n");
  for (int i = 0; i < 1; i++) // only a subset of points
    fprintf(stderr, "[x y z] = %.3lf %.3lf %.3lf s = %.3lf\n", 
	    x[i], y[i], z[i], sc[*scalar_id][i]);
  fprintf(stderr, "vtk cells\n");
  for (int i = 0; i < 1; i++) { // only print a subset of blocks
    int s = i * *nx * *ny * *nz;
    for (int j = 0; j < *nx * *ny * *nz; j++) {
      if ((j % *nx == *nx - 1) || ((j / *nx) % *ny == *ny - 1) ||
	  (*nz > 1 && (j / (*nx * *ny)) % *nz == *nz - 1))
	continue;
      fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j], 
	      y[s + j], z[s + j]);
      fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + 1],
	      y[s + j + 1], z[s + j + 1]);
      fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx + 1],
	      y[s + j + *nx + 1], z[s + j + *nx + 1]);
      fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx],
	      y[s + j + *nx], z[s + j + *nx]);
      if (*nz > 1) {
	fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx * *ny],
		y[s + j + *nx * *ny], z[s + j + *nx * *ny]);
	fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx * *ny + 1],
		y[s + j + *nx * *ny + 1], z[s + j + *nx * *ny + 1]);
	fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx + *nx * *ny + 1],
		y[s + j + *nx + *nx * *ny + 1], z[s + j + *nx + *nx * *ny + 1]);
	fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx + *nx * *ny],
		y[s + j + *nx + *nx * *ny], z[s + j + *nx + *nx * *ny]);
      }
      fprintf(stderr, "\n");
    }
  } // first two local blocks

  // init scalar field
  vtkDoubleArray *scalars = vtkDoubleArray::New();
  scalars->SetNumberOfValues(nv);
  double scalar_min = sc[*scalar_id][0];
  double scalar_max = sc[*scalar_id][0];
  for(int i = 0; i < nv; i++) {
    scalars->SetValue(i, sc[*scalar_id][i]);
    if (sc[*scalar_id][i] < scalar_min)
      scalar_min = sc[*scalar_id][i];
    if (sc[*scalar_id][i] > scalar_max)
      scalar_max = sc[*scalar_id][i];
    // debug
//     if (fabs(sc[*scalar_id][i]) > 2)
//       fprintf(stderr, "scalar[%d] = %.3e\n", i, sc[*scalar_id][i]);
  }
  // debug
  fprintf(stderr, "scalar range = %.3lf %.3lf\n", scalar_min, scalar_max);

  vtkUnstructuredGrid *data = vtkUnstructuredGrid::New();
  data->SetPoints(points);
  if (*nz == 1) // 2D
    data->SetCells(VTK_QUAD, cells);
  else // 3D
    data->SetCells(VTK_HEXAHEDRON, cells);
  data->GetPointData()->SetScalars(scalars);

  vtkDataSetMapper *mapper = vtkDataSetMapper::New();
  mapper->SetInput(data);
  mapper->SetScalarRange(scalar_min, scalar_max);

  // init actor
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);

  // init renderer
  vtkRenderer *renderer = vtkRenderer::New();
  renderer->AddActor(actor);
  renderer->SetBackground(0.1, 0.2, 0.4);

  // init window
  vtkRenderWindow *window = vtkRenderWindow::New();
  window->AddRenderer(renderer);
  window->SetSize(512, 512);

  // init window interaction and run the window
  vtkRenderWindowInteractor *interactor = vtkRenderWindowInteractor::New();
  if (rank == 0) { // only draw from one rank
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    interactor->Start();
  };

  // cleanup
  coords->Delete();
  points->Delete();
  cells->Delete();
  scalars->Delete();
  data->Delete();
  mapper->Delete();
  actor->Delete();
  renderer->Delete();
  window->Delete();
  interactor->Delete();

}
//--------------------------------------------------------------------------
