/*
 * ITL_geodesictree.h
 *
 *  Created on: Apr 18, 2012
 *      Author: abon
 */

/*
#ifndef ITL_GEODESICTREE_H_
#define ITL_GEODESICTREE_H_

#include "ITL_header.h"
#include "ITL_trianglepatch.h"

class ITL_geodesictree
{
	double PHI;											// Golden ratio.
	double NUM_TRIANGLES_IN_AN_ICOSAHEDRON;				// The number of faces in an icosahedron.
	int NUM_DIMENSIONS;									// The dimensionality of the sphere.
	int NUM_VERTICES_IN_A_TRIANGLE;						// The number of vertices in a triangle.
	int NUM_NEW_TRIANGLES_PER_TRIANGLE;					// The number of triangles produced

	list<VECTOR3> vertexList;
	//list<ITL_trianglepatch> trianglePatchList;

	ITL_trianglepatch* root;

public:

	ITL_geodesictree()
	{
		double PHI = 1.6180339887;
		double NUM_TRIANGLES_IN_AN_ICOSAHEDRON = 20;
		NUM_DIMENSIONS = 3;
		NUM_VERTICES_IN_A_TRIANGLE = 3;
		NUM_NEW_TRIANGLES_PER_TRIANGLE = 4;

	}

	void
	initTree()
	{
		// Add origin to vertexlist
		VECTOR3 origin( 0, 0, 0 );
		vertexList.push_back( origin );

		// Create root node (at origin)
		root = new ITL_trianglepatch( 0, 0, 0, 0 );

		double vertices[11][3];
		vertices[0][0] = 0; vertices[0][1] = PHI; vertices[0][2] = 1;
		vertices[1][0] = 0; vertices[1][1] = -PHI; vertices[1][2] = 1;
		vertices[2][0] = 0; vertices[2][1] = PHI; vertices[2][2] = -1;
		vertices[3][0] = 0; vertices[3][1] = -PHI; vertices[3][2] = -1;
		vertices[4][0] = 1; vertices[4][1] = 0; vertices[4][2] = PHI;
		vertices[5][0] = -1; vertices[5][1] = 0; vertices[5][2] = PHI;
		vertices[6][0] = 1; vertices[6][1] = 0; vertices[6][2] = -PHI;
		vertices[7][0] = PHI; vertices[7][1] = 1; vertices[7][2] = 0;
		vertices[8][0] = -PHI; vertices[8][1] = 1; vertices[8][2] = 0;
		vertices[9][0] = PHI; vertices[9][1] = -1; vertices[9][2] = 0;
		vertices[10][0] = -PHI; vertices[10][1] = -1; vertices[10][2] = 0;
		double mg = sqrt( 1*1 + PHI*PHI );
		for( int i=0; i<11; i++ )
		{
			vertices[i][0] /= mg;
			vertices[i][1] /= mg;
			vertices[i][2] /= mg;
			VECTOR3 nextVertex( vertices[i][0], vertices[i][1], vertices[i][2] );
			vertexList.push_back( nextVertex );

		}

		double A[20] = { 1, 4, 8, 6, 10, 3, 5, 1, 0, 4, 2, 8, 7, 6, 11, 11, 5, 0, 2, 9 };
		double B[20] = { 3, 10, 4, 8, 6, 1, 11, 4, 5, 8, 0, 6, 2, 3, 7, 5, 0, 2, 9, 7 };
		double C[20] = { 10, 1, 10, 10, 3, 11, 1, 5, 4, 0, 8, 2, 6, 7, 3, 9, 9, 9, 7, 11 };

		for( int i=0; i<20; i++ )
		{
			ITL_trianglepatch nextpatch( A[i], B[i], C[i], 0 );
			trianglePatchList.push_back( nextpatch );
		}

	}

	void
	createTree( )
	{

	}

};
*/


//#endif /* ITL_GEODESICTREE_H_ */
