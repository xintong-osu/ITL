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
#include "ITL_util.h"

#include <vector>

#include "iMesh.h"
#include "MBiMesh.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "moab/ParallelComm.hpp"

// vtk headers
#include "vtkProperty.h"
#include  "vtkLookupTable.h"
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
#include "vtkArrowSource.h"
#include "vtkGlyph2D.h"
#include "vtkGlyph3D.h"
#include "vtkSmartPointer.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"

using namespace std;
using namespace moab;

// vortex
struct vortex_t {
  double x, y, z; 	// location
  int t; 		// time
  double l2; 		// lambda-2 value
};

static int timestep = 0;
static int nb; 			// number of local blocks
static int nv; 			// number of local vertices in all local blocks
static int rank = 0; 		// MPI rank
static int groupsize = 1; 	// MPI communicator size

static double *x, *y, *z; 	// Vertex geometry
static double **sc;			// Vertex property of scalar fields
static double **press;		// Pressure of scalar fields
static double **vx, **vy, **vz; // Vector fields, with separate components
static int nsf, nvf; 		// number of scalar and vector fields

static vector<vortex_t> vortices;

static int nElement = 0;
static int nVertexElement = 0;
static int selectedElementId = 0;
vtkSmartPointer<vtkDoubleArray> coords;
vtkSmartPointer<vtkPoints> points;
vtkSmartPointer<vtkDoubleArray> scalars;
vtkSmartPointer<vtkCellArray> cells;
vtkSmartPointer<vtkUnstructuredGrid> data;
vtkSmartPointer<vtkDataSetMapper> mapper;
vtkSmartPointer<vtkLookupTable> colorLookupTable;
vtkSmartPointer<vtkUnsignedCharArray> colors;
vtkSmartPointer<vtkActor> actor;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> window;
vtkSmartPointer<vtkRenderWindowInteractor> interactor;
vtkSmartPointer<vtkCallbackCommand> keypressCallback;


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

  press = new double*[nsf];
  sc = new double*[nsf];
  for (int i = 0; i < nsf; i++)
    press[i] = new double[nv];
  for (int i = 0; i < nsf; i++)
    sc[i] = new double[nv];

  vx = new double*[nvf];
  vy = new double*[nvf];
  vz = new double*[nvf];
  for (int i = 0; i < nvf; i++) {
    vx[i] = new double[nv];
    vy[i] = new double[nv];
    vz[i] = new double[nv];
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &groupsize);

  if (rank == 0)
    fprintf( stderr, "Beginning ITL. Number of blocks = %d "
    				 "Number of vertices per block = %d "
    				 "Number of scalar fields = %d "
    				 "Number of vector fields = %d\n",
    				 nb, *num_verts_per_block, nsf, nvf );

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
    fprintf(stderr, "Ending ITL\n");

}

//--------------------------------------------------------------------------
// 
// get geometry
//
// x, y, z: arrays of x coordinates, y coordinates, z coordinates for
//   cell positions (meaning of a cell position, eg. cell center, vertex, etc.
//   depends on the application. CIAN considers the scalar values to be
//   at the provided coordinates
// npts: number of x, y, z points
// ofst: starting offset (index) of points
//

// C and C++ wrapper
void
ITL_geom()
{
}

// fortran wrapper
extern "C"
void
itl_geom_( double *xx, double *yy, double *zz, int *npts, int *ofst )
{
	assert( x != NULL );
	assert( y != NULL );
	assert( z != NULL );
	assert( xx != NULL );
	assert( yy != NULL );
	assert( zz != NULL );
	assert( npts != NULL );
	assert( ofst != NULL );

	fprintf( stderr, "ITL: geometry %d points\n", *npts );
	fprintf( stderr, "ITL: offset %d\n", *ofst );

	for (int i = 0; i < *npts; i++)
	{
		fprintf( stderr, "%d\n", i );
		x[*ofst + i] = xx[*ofst + i];
		y[*ofst + i] = yy[*ofst + i];
		z[*ofst + i] = zz[*ofst + i];
	}

	fprintf( stderr, "ITL: geometry %d points\n", *npts );

	// debug
	if (rank == 0)
	{
		fprintf(stderr, "ITL: geometry %d points:\n", *npts);
		for (int i = 0; i < 64; i++)
			fprintf(stderr, "[x y z] = %.3lf %.3lf %.3lf\n",
							x[*ofst + i], y[*ofst + i], z[*ofst + i]);
	}

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
void ITL_scalar_data()
{
}

// fortran wrapper
extern "C"
void
itl_scalar_data_(int *scalar_id, double *vals, int *nvals, int *ofst)
{
	fprintf( stderr, "In itl_scalar_data\n" );
	fprintf( stderr, "Number of vertices: %d %d\n", *nvals, *ofst );

	// Set number of blocks / fine grain elements
	nElement = (*nvals);

	// Copy data
	for (int i = 0; i < *nvals; i++)
		for (int j = 0; j < *ofst; j++)
			sc[*scalar_id][i * (*ofst) + j] = vals[i * (*ofst) + j];
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

    fprintf(stderr, "ITL: vector field:\n");
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

//--------------------------------------------------------------------------
void
keyPressCallbackFunc( vtkObject* caller,
					  long unsigned int vtkNotUsed(eventId),
					  void* vtkNotUsed(clientData),
					  void* vtkNotUsed(callData) )
{
	double dcolor[3];
	unsigned char color[3];

	fprintf( stderr, "Keypress detected: %s..\n", interactor->GetKeySym() );
	fprintf( stderr, "Selected element ID: %d..\n", selectedElementId );
	fprintf( stderr, "1\n" );

	// Put color of the previously selectes element back
	int startId = selectedElementId * nVertexElement;
	fprintf( stderr, "2\n" );
	for(int i = startId; i < (startId + nVertexElement); i++ )
	{
		// Get pressure at current point
		// and use that to decide coolr
		//colorLookupTable->GetColor( sc[0][i], dcolor );
		colorLookupTable->GetColor( selectedElementId, dcolor );
		//p = sqrt( vx[0][i]*vx[0][i] + vy[0][i]*vy[0][i] + vz[0][i]*vz[0][i] );
		//colorLookupTable->GetColor( p, dcolor );

		#ifdef DEBUG_MODE
		fprintf( stderr, "dcolor: <%g %g %g>\n", dcolor[0], dcolor[1], dcolor[2] );
		#endif

		for(unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}

		#ifdef DEBUG_MODE
		fprintf( stderr, "color: <%d %d %d>\n", (int)color[0], (int)color[1], (int)color[2] );
		#endif

		colors->SetTupleValue( i, color );
	}
	fprintf( stderr, "3\n" );

	// Increment selected block ID
	fprintf( stderr, "Selected element ID: %d..\n", selectedElementId );
	selectedElementId  = (selectedElementId + 1) % nElement;
	fprintf( stderr, "Selected element ID: %d..\n", selectedElementId );
	fprintf( stderr, "4\n" );

	// Change color of the currently selected element to white
	startId = selectedElementId * nVertexElement;
	for(int i = startId; i < (startId + nVertexElement); i++ )
	{
		for(unsigned int j = 0; j < 3; j++)
			color[j] = static_cast<unsigned char>( 255.0 );

		#ifdef DEBUG_MODE
		fprintf( stderr, "color: <%d %d %d>\n", (int)color[0], (int)color[1], (int)color[2] );
		#endif

		colors->SetTupleValue( i, color );
	}
	fprintf( stderr, "5\n" );

	window->Render();
}

//--------------------------------------------------------------------------
//
// renders a scalar field using vtk
//
// scalar_id: index of scalar field [0 - num_scalar_fields - 1]
// nx, ny, nz: block size
//
// C and C++ wrapper
void
ITL_vis_scalar()
{
}

// fortran wrapper
extern "C"
void
itl_vis_scalar_(int *scalar_id, int *nx, int *ny, int *nz)
{
	// init vtk float array
	coords = vtkSmartPointer<vtkDoubleArray>::New();
	coords->SetNumberOfComponents(3);
	coords->SetNumberOfTuples(nv);
	for (int i = 0; i < nv; i++)
	{
		//if( i < 100 )
		//	fprintf( stderr, "Test: %g %g %g\n", x[i], y[i], z[i] );
		coords->SetTuple3(i, x[i], y[i], z[i]);
	}

	// init vtk points
	points = vtkSmartPointer<vtkPoints>::New();
	points->SetData(coords);

	// init vtk cells
	cells = vtkSmartPointer<vtkCellArray>::New();
	for (int i = 330; i < 340; i++)
	{
		// for all local blocks
		int s = i * *nx * *ny * *nz;

		for (int j = 0; j < *nx * *ny * *nz; j++)
		{
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

			if (*nz > 1)
			{
				cells->InsertCellPoint(s + j + *nx * *ny);
				cells->InsertCellPoint(s + j + *nx * *ny + 1);
				cells->InsertCellPoint(s + j + *nx + *nx * *ny + 1);
				cells->InsertCellPoint(s + j + *nx + *nx * *ny);
			}

		} // inner for
	} // outer for : local blocks

	// Debug
	#ifdef DEBUG_MODE
	fprintf( stderr, "vtk scalar field\n" );
	for (int i = 0; i < 1; i++) // only a subset of points
		fprintf(stderr, "[x y z] = %.3lf %.3lf %.3lf s = %.3lf\n",
						x[i], y[i], z[i], sc[*scalar_id][i] );

	fprintf(stderr, "vtk cells\n");
	for (int i = 0; i < 1; i++)
	{
		// only print a subset of blocks
		int s = i * *nx * *ny * *nz;
		for (int j = 0; j < *nx * *ny * *nz; j++)
		{
			if ((j % *nx == *nx - 1) || ((j / *nx) % *ny == *ny - 1) ||
				(*nz > 1 && (j / (*nx * *ny)) % *nz == *nz - 1))
				continue;
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j], y[s + j], z[s + j]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + 1], y[s + j + 1], z[s + j + 1]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx + 1], y[s + j + *nx + 1], z[s + j + *nx + 1]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx], y[s + j + *nx], z[s + j + *nx]);

			if (*nz > 1)
			{
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
	#endif

	// init scalar field
	scalars = vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(nv);
	double scalar_min = sc[*scalar_id][0];
	double scalar_max = sc[*scalar_id][0];
	for(int i = 0; i < nv; i++)
	{
		scalars->SetValue( i, sc[*scalar_id][i] );
		if( sc[*scalar_id][i] < scalar_min )
			scalar_min = sc[*scalar_id][i];
		if( sc[*scalar_id][i] > scalar_max )
			scalar_max = sc[*scalar_id][i];

		#ifdef DEBUG_MODE
		     if (fabs(sc[*scalar_id][i]) > 2)
		       fprintf(stderr, "scalar[%d] = %.3e\n", i, sc[*scalar_id][i]);
		#endif
	}

	// Debug
	fprintf( stderr, "Number of elements = %d\n", nElement );
	fprintf( stderr, "Number of vertices/elements = %d\n", nVertexElement );
	fprintf( stderr, "Number of vertices = %d\n", nv );
	fprintf( stderr, "Scalar range = %g %g\n", scalar_min, scalar_max );

	data = vtkSmartPointer<vtkUnstructuredGrid>::New();
	data->SetPoints( points );
	if ( *nz == 1 ) // 2D
		data->SetCells(VTK_QUAD, cells);
	else // 3D
		data->SetCells(VTK_HEXAHEDRON, cells);
	data->GetPointData()->SetScalars( scalars );

	// Create the color map
	colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
	//
	//colorLookupTable->SetTableRange( scalar_min, scalar_max );
	// OR
	colorLookupTable->SetTableRange( 0.0, 1.0 );
	//
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	double p;
	double dcolor[3];
	unsigned char color[3];

	// Color vertices by pressure value
	/*
	for(int i = 0; i <nv; i++ )
	{
		// Get pressure at current point
		// and use that to decide coolr
		colorLookupTable->GetColor( sc[*scalar_id][i], dcolor );

		//p = sqrt( vx[0][i]*vx[0][i] + vy[0][i]*vy[0][i] + vz[0][i]*vz[0][i] );
		//colorLookupTable->GetColor( p, dcolor );

		#ifdef DEBUG_MODE
		fprintf( stderr, "dcolor: <%g %g %g>\n", dcolor[0], dcolor[1], dcolor[2] );
		#endif

		for(unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}

		#ifdef DEBUG_MODE
		fprintf( stderr, "color: <%d %d %d>\n", (int)color[0], (int)dcolor[1], (int)dcolor[2] );
		#endif

		colors->InsertNextTupleValue( color );
	}
	*/

	// Color blocks by block ID
	// OR alternate between two colors
	for(int i = 0; i <nElement; i++ )
	{
		for( int j = 0; j<nVertexElement; j++ )
		{
			//colorLookupTable->GetColor( i, dcolor );
			// OR
			if( i%2 == 0 )
				colorLookupTable->GetColor( 0.2, dcolor );
			else
				colorLookupTable->GetColor( 0.8, dcolor );

			#ifdef DEBUG_MODE
			fprintf( stderr, "dcolor: <%g %g %g>\n", dcolor[0], dcolor[1], dcolor[2] );
			#endif

			for(unsigned int k = 0; k < 3; k++)
			{
				color[k] = static_cast<unsigned char>(255.0 * dcolor[k]);
			}

			#ifdef DEBUG_MODE
			fprintf( stderr, "color: <%d %d %d>\n", (int)color[0], (int)dcolor[1], (int)dcolor[2] );
			#endif

			colors->InsertNextTupleValue( color );
		}
	}

	// Color blocks randomly
	/*
	srand ( time(NULL) );
	for(int i = 0; i <nElement; i++ )
	{
		for( int j = 0; j<nVertexElement; j++ )
		{
			// Get pressure at current point
			// and use that to decide coolr
			colorLookupTable->GetColor( ( rand() % 100 ) / 99.0, dcolor );

			#ifdef DEBUG_MODE
			fprintf( stderr, "dcolor: <%g %g %g>\n", dcolor[0], dcolor[1], dcolor[2] );
			#endif

			for(unsigned int k = 0; k < 3; k++)
			{
				color[k] = static_cast<unsigned char>( 255.0 * dcolor[k] );
			}

			#ifdef DEBUG_MODE
			fprintf( stderr, "color: <%d %d %d>\n", (int)color[0], (int)dcolor[1], (int)dcolor[2] );
			#endif

			colors->InsertNextTupleValue( color );
		}
	}
	*/

	// Load the colors to the vertices
	data->GetPointData()->SetScalars(colors);

	mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInput(data);
	//mapper->SetInputData( data );
	mapper->SetScalarRange(scalar_min, scalar_max);

	// init actor
	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToWireframe();

	// init renderer
	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(0.1, 0.2, 0.4);

	// init window
	window = vtkSmartPointer<vtkRenderWindow>::New();
	window->AddRenderer(renderer);
	window->SetSize(512, 512);

	// init window interaction and run the window
	interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	// init keyboard callback
	keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	keypressCallback->SetCallback ( keyPressCallbackFunc );
	interactor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

	// Start rendering
	if (rank == 0)
	{
		// only draw from one rank
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

//--------------------------------------------------------------------------
//
// renders a vector field using vtk
//
// scalar_id: index of scalar field [0 - num_scalar_fields - 1]
// nx, ny, nz: block size
//
// C and C++ wrapper
void
ITL_vis_vector()
{
}

// fortran wrapper
extern "C"
void
itl_vis_vector_(int *scalar_id, int *nx, int *ny, int *nz)
{
	// Init vtk float array
	coords = vtkSmartPointer<vtkDoubleArray>::New();
	coords->SetNumberOfComponents(3);
	coords->SetNumberOfTuples(nv);
	for (int i = 0; i < nv; i++)
		coords->SetTuple3(i, x[i], y[i], z[i]);

	// Init vtk points
	points = vtkSmartPointer<vtkPoints>::New();
	points->SetData(coords);

	// init vtk cells
	cells = vtkSmartPointer<vtkCellArray>::New();
	for (int i = 0; i < nb; i++)
	{
		// for all local blocks
		int s = i * *nx * *ny * *nz;

		for (int j = 0; j < *nx * *ny * *nz; j++)
		{
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

			if (*nz > 1)
			{
				cells->InsertCellPoint(s + j + *nx * *ny);
				cells->InsertCellPoint(s + j + *nx * *ny + 1);
				cells->InsertCellPoint(s + j + *nx + *nx * *ny + 1);
				cells->InsertCellPoint(s + j + *nx + *nx * *ny);
			}

		} // inner for
	} // outer for : local blocks

	// Debug
	#ifdef DEBUG_MODE
	fprintf( stderr, "vtk vector field\n" );
	for (int i = 0; i < 1; i++) // only a subset of points
		fprintf(stderr, "[x y z] = %.3lf %.3lf %.3lf\n",
						x[i], y[i], z[i] );

	fprintf(stderr, "vtk cells\n");
	for (int i = 0; i < 1; i++)
	{
		// only print a subset of blocks
		int s = i * *nx * *ny * *nz;
		for (int j = 0; j < *nx * *ny * *nz; j++)
		{
			if ((j % *nx == *nx - 1) || ((j / *nx) % *ny == *ny - 1) ||
				(*nz > 1 && (j / (*nx * *ny)) % *nz == *nz - 1))
				continue;
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j], y[s + j], z[s + j]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + 1], y[s + j + 1], z[s + j + 1]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx + 1], y[s + j + *nx + 1], z[s + j + *nx + 1]);
			fprintf(stderr, "[%.3lf %.3lf %.3lf] ", x[s + j + *nx], y[s + j + *nx], z[s + j + *nx]);

			if (*nz > 1)
			{
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
	#endif

	// Init scalar field
	scalars = vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfValues(nv);
	double scalar_min = sc[*scalar_id][0];
	double scalar_max = sc[*scalar_id][0];
	for(int i = 0; i < nv; i++)
	{
		scalars->SetValue( i, sc[*scalar_id][i] );
		if( sc[*scalar_id][i] < scalar_min )
			scalar_min = sc[*scalar_id][i];
		if( sc[*scalar_id][i] > scalar_max )
			scalar_max = sc[*scalar_id][i];
		// debug
		//     if (fabs(sc[*scalar_id][i]) > 2)
		//       fprintf(stderr, "scalar[%d] = %.3e\n", i, sc[*scalar_id][i]);
	}

	// Setup the arrows
	vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
	arrowSource->Update();

	vtkSmartPointer<vtkGlyph3D> glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
	glyphFilter->SetSourceConnection( arrowSource->GetOutputPort() );
	glyphFilter->OrientOn();
	glyphFilter->SetVectorModeToUseVector();

	//glyphFilter->SetInputData( coords );
	glyphFilter->Update();

	// Debug
	fprintf(stderr, "scalar range = %.3lf %.3lf\n", scalar_min, scalar_max);

	data = vtkSmartPointer<vtkUnstructuredGrid>::New();
	data->SetPoints(points);
	if (*nz == 1) // 2D
		data->SetCells(VTK_QUAD, cells);
	else // 3D
		data->SetCells(VTK_HEXAHEDRON, cells);
	data->GetPointData()->SetScalars(scalars);

	mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInput(data);
	mapper->SetScalarRange(scalar_min, scalar_max);

	// init actor
	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// init renderer
	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(0.1, 0.2, 0.4);

	// init window
	window = vtkSmartPointer<vtkRenderWindow>::New();
	window->AddRenderer(renderer);
	window->SetSize(512, 512);

	// init window interaction and run the window
	interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	if (rank == 0)
	{
		// only draw from one rank
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
int
array_index_3dto1d( int x, int y, int z, int* dim )
{
	if( x < 0 ) x = 0;
	if( y < 0 ) y = 0;
	if( z < 0 ) z = 0;
	if( x > dim[0]-1 ) x = dim[0]-1;
	if( y > dim[1]-1 ) y = dim[1]-1;
	if( z > dim[2]-1 ) z = dim[2]-1;

	return ( z * dim[0]*dim[1] + y * dim[0] + x );
}

//--------------------------------------------------------------------------

// C and C++ wrapper
void
ITL_resample()
{
}

// fortran wrapper
extern "C"
void
itl_resample_( int count, int tot_pts )
{
	double xmin = 100000, xmax = -100000;
	double ymin = 100000, ymax = -100000;
	double zmin = 100000, zmax = -100000;

	for( int i = 0; i < count; i++ )
	{
		for (int j = 0; j < tot_pts; j++)
		{
			if( x[i * tot_pts + j] < xmin ) xmin = x[i * tot_pts + j];
			if( x[i * tot_pts + j] > xmax ) xmax = x[i * tot_pts + j];
			if( y[i * tot_pts + j] < ymin ) ymin = y[i * tot_pts + j];
			if( y[i * tot_pts + j] > ymax ) ymax = y[i * tot_pts + j];
			if( z[i * tot_pts + j] < zmin ) zmin = z[i * tot_pts + j];
			if( z[i * tot_pts + j] > zmax ) zmax = z[i * tot_pts + j];

		}
	}
	fprintf( stderr, "Range: <%g %g, %g %g, %g %g>\n", xmin, xmax, ymin, ymax, zmin, zmax );

	// Determine spatial extent of the field
	int fieldDim[3];
	double xRange = xmax - xmin;
	double yRange = ymax - ymin;
	double zRange = zmax - zmin;
	fieldDim[2] = 128;
	fieldDim[0] = (int) floor( ( xRange / zRange ) * fieldDim[2] );
	fieldDim[1] = (int) floor( ( yRange / zRange ) * fieldDim[2] );
	fprintf( stderr, "Field Size: <%d %d %d>\n", fieldDim[0], fieldDim[1], fieldDim[2] );

	// Allocate memory
	int nTotPoint = fieldDim[0] * fieldDim[1] * fieldDim[2];
	float* resampledField = new float[nTotPoint];
	memset( resampledField, 0, sizeof(float)*nTotPoint );

	// Compute pressure at each point
	double xloc, yloc, zloc, pr, velx, vely, velz, velmag;
	double sumScalar = 0;
	double dist = 0;
	double w = 0;
	int ind = 0;
	double minDist = 1000;
	double curx, cury, curz;
	double deltax = xRange / (fieldDim[0] - 1);
	double deltay = yRange / (fieldDim[1] - 1);
	double deltaz = zRange / (fieldDim[2] - 1);

	double x_nbh = deltax;
	double y_nbh = deltay;
	double z_nbh = deltaz;
	double decay;
	float m, M;

	for( int r = 0; r<fieldDim[2]; r++ )
	{
		// Set current z-location
		curz = zmin + r * deltaz;

		for( int q = 0; q<fieldDim[1]; q++ )
		{
			// Set current y-location
			cury = ymin + q * deltay;

			for( int p = 0; p<fieldDim[0]; p++ )
			{
				// Set current x-location
				curx = xmin + p * deltax;

				sumScalar = 0;
				w = 0;
				minDist = 1000;

				// Find nearby irregular grid points and average their pressure
				for( int iC = 0; iC < count; iC++ )
				{
					for (int jP = 0; jP < tot_pts; jP++)
					{
						xloc = x[iC * tot_pts + jP];
						yloc = y[iC * tot_pts + jP];
						zloc = z[iC * tot_pts + jP];

						//pr = press[0][iC * tot_pts + jP];

						velx = vx[0][iC * tot_pts + jP];
						vely = vy[0][iC * tot_pts + jP];
						velz = vz[0][iC * tot_pts + jP];
						velmag = sqrt( velx*velx + vely*vely + velz*velz );

						if( (curx - xloc) <= x_nbh &&
							(cury - yloc) <= y_nbh &&
							(curz - zloc) <= z_nbh
						 )
						//if( dist < 1.0  )
						{
							dist = sqrt( (curx - xloc)*(curx - xloc) +
										 (cury - yloc)*(cury - yloc) +
										 (curz - zloc)*(curz - zloc) );

							//decay = exp( -dist );
							if( dist == 0 )
								decay = 1;
							else
								decay = 1.0 / (dist*dist);

							//sumScalar += pr * decay;
							sumScalar += velmag * decay;

							w += decay;

							if( dist == 0 )
							{
								fprintf( stderr, "%g %g %g %g %g %g\n", curx, cury, curz, xloc, yloc, zloc );
								fprintf( stderr, "%g %g %g %g\n", dist, decay, sumScalar, w );
							}


						}// end if
					}
				}

				// Get average pressure
				if( w == 0 )
					resampledField[ind] = 0;
				else
					resampledField[ind] = (float) ( sumScalar / w );

				#ifdef DEBUG_MODE
				if( ind == 0 )
				{
					m = resampledField[ind];
					M = resampledField[ind];
				}
				else
				{
					if( resampledField[ind] < m ) m = resampledField[ind];
					if( resampledField[ind] > M ) M = resampledField[ind];
				}
				#endif

				// Move to next grid point
				ind ++;

			}// end p : grid x
		}// end q : grid y
	}// end r : grid z

	#ifdef DEBUG_MODE
	fprintf( stderr, "Scalar field min/max: %g %g\n", m, M );
	#endif

	// Save pressure Field
	#ifdef DEBUG_MODE
	FILE* pressurefile = fopen( "/homes/chaudhua/nek_pressure.vol", "wb" );
	fwrite( fieldDim, sizeof(int), 3, pressurefile );
	fwrite( resampledField, sizeof(float)*nTotPoint, 1, pressurefile );
	fclose( pressurefile );

	FILE* pressurefile2 = fopen( "/homes/chaudhua/nek_pressure.csv", "w" );
	for( int i=0; i<nTotPoint; i++ )
		fprintf( pressurefile2, "%g\n", resampledField[i] );
	fclose( pressurefile2 );
	#endif

	// Save velocity magnitude field
	#ifdef DEBUG_MODE
	FILE* velmagfile = fopen( "/homes/chaudhua/nek_velmag.vol", "wb" );
	fwrite( fieldDim, sizeof(int), 3, velmagfile );
	fwrite( resampledField, sizeof(float)*nTotPoint, 1, velmagfile );
	fclose( velmagfile );
	#endif

	// Free memory
	delete [] resampledField;

}

// C and C++ wrapper
void
ITL_mesh()
{
}

// fortran wrapper
extern "C"
void
itl_mesh_(iMesh_Instance *mesh)
{
	fprintf( stderr, "Current time: %d\n", timestep );

	if( timestep < 50 )
	{
		timestep ++;
		return;
	}

	int ierr;
	Tag tag_x, tag_y, tag_z, tag_dims;
	Tag tag_vx, tag_vy, tag_vz, tag_press;
	int dims[3];
	ErrorCode merr;
	Range sem_sets;
	Range blocks;

	fprintf(stderr, "in function ITL_mesh ...\n");

	// Only printing rank 0 to reduce the amount of output
	if( rank == 0 )
	{
		for (int dim = iBase_VERTEX; dim <= iBase_REGION; dim++)
		{
			int numd;
			iMesh_getNumOfType(*mesh, 0, dim, &numd, &ierr);
			fprintf(stderr, "Number of %dd elements = %d\n", dim, numd);
		}
	}

	// cast the imesh instance to a moab instance
	moab::Interface *mb = ((MBiMesh *)*mesh)->mbImpl;

	EntityHandle root_set = mb->get_root_set();

	// Get block size along each dimension
	merr = mb->tag_get_handle( "SEM_DIMS", 3, MB_TYPE_INTEGER, tag_dims );
	fprintf( stderr, "%d\n", merr );
	assert(merr == moab::MB_SUCCESS );

	merr = mb->get_entities_by_type_and_tag( root_set, moab::MBENTITYSET, &tag_dims, NULL, 1, sem_sets );
	assert(merr == MB_SUCCESS);

	EntityHandle sem_set = sem_sets[0];
	merr = mb->tag_get_data( tag_dims, &sem_set, 1, (void*)dims );
	assert(merr == MB_SUCCESS);
	fprintf(stderr, "dims = [%d %d %d]\n", dims[0], dims[1], dims[2]);

	// get blocks
	merr = mb->get_entities_by_type(0, MBHEX, blocks);
	assert(merr == MB_SUCCESS);

	//merr = mb->tag_get_data( tag_dims, &root_set, 1, dims );
	//fprintf( stderr, "%d\n", merr );
	//fprintf( stderr, "%d %d %d\n", dims[0], dims[1], dims[2] );
	//assert(merr == moab::MB_SUCCESS );
	dims[0] = dims[1] = dims[2] = 4;

	// get the tag handles for the x, y, z geometry arrays
	nVertexElement = dims[0] * dims[1] * dims[2];

	merr = mb->tag_get_handle("SEM_X", nVertexElement, MB_TYPE_DOUBLE, tag_x);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("SEM_Y", nVertexElement, MB_TYPE_DOUBLE, tag_y);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("SEM_Z", nVertexElement, MB_TYPE_DOUBLE, tag_z);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("PRESS", nVertexElement, MB_TYPE_DOUBLE, tag_press);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("VX", nVertexElement, MB_TYPE_DOUBLE, tag_vx);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("VY", nVertexElement, MB_TYPE_DOUBLE, tag_vy);
	assert(merr == MB_SUCCESS);

	merr = mb->tag_get_handle("VZ", nVertexElement, MB_TYPE_DOUBLE, tag_vz);
	assert(merr == MB_SUCCESS);


	// get x,y,z vertex geometry for each block
	int count;

	mb->tag_iterate( tag_x, blocks.begin(), blocks.end(), count, (void *&)x );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "x-values read\n" );

	mb->tag_iterate( tag_y, blocks.begin(), blocks.end(), count, (void *&)y );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "y-values read\n" );

	mb->tag_iterate( tag_z, blocks.begin(), blocks.end(), count, (void *&)z );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "z-values read\n" );

	mb->tag_iterate(tag_press, blocks.begin(), blocks.end(), count, (void *&)press[0] );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "pressure values read\n" );

	mb->tag_iterate(tag_vx, blocks.begin(), blocks.end(), count, (void *&)vx[0] );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "vx values read\n" );

	mb->tag_iterate(tag_vy, blocks.begin(), blocks.end(), count, (void *&)vy[0] );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "vy values read\n" );

	mb->tag_iterate(tag_vz, blocks.begin(), blocks.end(), count, (void *&)vz[0] );
	assert(merr == MB_SUCCESS);
	assert(count == (int)blocks.size());

	fprintf(stderr, "vz values read\n" );

	// Print and write out vertex geometry
	#ifdef DEBUG_MODE
	fprintf(stderr, "ITL: MOAB geometry %d blocks\n", count);
	FILE* meshfile = fopen( "/homes/chaudhua/nek_mesh.csv", "w" );
	double p_min = press[0][0];
	double p_max = press[0][0];
	for( int i = 0; i < count; i++ )
	{
	//for (int i = 0; i < 1; i++) { // print only 1 block for now
		for (int j = 0; j < tot_pts; j++)
		{
			//fprintf( stderr, "[x y z] = %.3lf, %.3lf, %.3lf "
			//				   "[vx vy vz] = %.3lf, %.3lf, %.3lf "
			//				   "[press] = %.3lf\n",
			//				   x[i * tot_pts + j], y[i * tot_pts + j],	z[i * tot_pts + j],
			//				   vx[i * tot_pts + j], vy[i * tot_pts + j],vz[i * tot_pts + j],
			//				   press[i * tot_pts + j] );
			fprintf( meshfile, "%g, %g, %g, %g, %g, %g, %g\n",
							   x[i * tot_pts + j], y[i * tot_pts + j],	z[i * tot_pts + j],
							   vx[0][i * tot_pts + j], vy[0][i * tot_pts + j],vz[0][i * tot_pts + j],
							   press[0][i * tot_pts + j] );

			//if( press[0][i * tot_pts + j] < p_min ) p_min = press[0][i * tot_pts + j];
			//if( press[0][i * tot_pts + j] > p_max ) p_max = press[0][i * tot_pts + j];
		}
	}
	fclose( meshfile );
	//fprintf( stderr, "Pressure min/max: %g %g\n", p_min, p_max );
	#endif

	#ifdef DEBUG_MODE
	FILE* meshfile2 = fopen( "/homes/chaudhua/nek_mesh.bin", "wb" );
	fwrite( x, sizeof(double), count*tot_pts, meshfile2 );
	fwrite( y, sizeof(double), count*tot_pts, meshfile2 );
	fwrite( z, sizeof(double), count*tot_pts, meshfile2 );
	fwrite( vx[0], sizeof(double), count*tot_pts, meshfile2 );
	fwrite( vy[0], sizeof(double), count*tot_pts, meshfile2 );
	fwrite( vz[0], sizeof(double), count*tot_pts, meshfile2 );
	fwrite( press[0], sizeof(double), count*tot_pts, meshfile2 );
	fclose( meshfile2 );
	#endif

	// Resample vertex geometry
	//fprintf(stderr, "ITL: entering itl_resample\n" );
	//itl_resample_( count, tot_pts );
	//fprintf(stderr, "ITL: exiting itl_resample\n" );

	// Copy data for rendering
	fprintf(stderr, "ITL: entering itl_scalar_data\n" );
	int ofst = 0;
	int sc_id = 0;
	itl_scalar_data_( &sc_id, press[0], &count, &nVertexElement );
	fprintf(stderr, "ITL: exiting itl_scalar_data\n" );

	// Render
	fprintf(stderr, "ITL: entering itl_vis_scalar\n" );
	itl_vis_scalar_( &sc_id, &dims[0], &dims[1], &dims[2] );
	fprintf(stderr, "ITL: exiting itl_vis_scalar\n" );

}
