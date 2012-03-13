/**
 * @file MainIT_irregularScalar.cpp
 * Application program for entropy computation of irregular scalar field.
 * Created on: Feb 23, 2012
 * @author Cong Wang
 */

#include <mpi.h>
#include "ITL_header.h"
#include "ITL_base.h"
#include "ITL_util.h"
#include "ITL_ioutil.h"
#include "ITL_vectormatrix.h"
#include "ITL_histogramconstants.h"
#include "ITL_localentropy.h"
#include "ITL_globalentropy.h"
#include "ITL_globaljointentropy.h"

#include "ITL_cell.h"
#include "ITL_field_unstructured.h"

#include <limits>

using namespace std;

int numProcs;
int myId;

int verboseMode = 1;
int nBin = 1000;
double execTime[5];
clock_t starttime, endtime;
ITL_histogram *histogram = NULL;
ITL_localentropy<SCALAR> *localEntropyComputer = NULL;
ITL_globalentropy<SCALAR> *globalEntropyComputer = NULL;

#define LOCAL_ENTROPY

void compute_localentropy_serial(ITL_grid_tetrahedral<SCALAR>* tetGrid)
{
	tetGrid->buildBoundingBoxInfor();
	tetGrid->buildTetAdjacentInfor();

	//vector<int>* intersectTets = new vector<int>[nVert];
	//vector<int>* containTets = new vector<int>[nVert];
	//tetGrid->getBoxIntersecTets(intersectTets, containTets);
	tetGrid->buildBoxIntersecTets();

	ITL_field_unstructured<SCALAR>* scalarField = new ITL_field_unstructured<SCALAR>(tetGrid);
	localEntropyComputer = new ITL_localentropy<SCALAR>( scalarField, histogram );
	localEntropyComputer->computeLocalEntropyOfField_Unstructured( nBin, false );
}

void compute_globalentropy_serial(ITL_grid_tetrahedral<SCALAR>* tetGrid)
{
	ITL_field_unstructured<SCALAR>* scalarField = new ITL_field_unstructured<SCALAR>(tetGrid);
	globalEntropyComputer = new ITL_globalentropy<SCALAR>( scalarField, histogram );
	globalEntropyComputer->computeGlobalEntropyOfField_Unstructured( nBin, false );
	printf("global entropy: %f\n", globalEntropyComputer->getGlobalEntropy());
}

int main( int argc, char** argv )
{
	// Initialize MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myId );

	const char* tetFile = "blunt.scalar";
	SCALAR* vlist = NULL;
	int* tlist = NULL;
	int nVert, nTet;
	SCALAR* scalarFieldData = ITL_ioutil<SCALAR>::readTetrahedralSerial( tetFile, vlist, tlist, nVert, nTet);

	ITL_grid_tetrahedral<SCALAR>* tetGrid = new ITL_grid_tetrahedral<SCALAR>(nVert, nTet);
	tetGrid->radius = 0.0217;
	
	//load vertexList
	for (int i = 0; i < nVert; ++i)
	{
		tetGrid->vertexList[i].x = vlist[4*i];
		tetGrid->vertexList[i].y = vlist[4*i + 1];
		tetGrid->vertexList[i].z = vlist[4*i + 2];
		tetGrid->vertexList[i].f = vlist[4*i + 3];
	}
	delete[] vlist;
	//load cellList
	for (int i = 0; i < nTet; ++i)
	{
		tetGrid->cellList[i].index = i;
		for (int j = 0; j < 4; ++j)
		{
			tetGrid->cellList[i].v[j] = tlist[4*i+j];
		}
	}
	delete[] tlist;

#ifndef LOCAL_ENTROPY
	compute_globalentropy_serial(tetGrid);
#else
	compute_localentropy_serial(tetGrid);
#endif
}

