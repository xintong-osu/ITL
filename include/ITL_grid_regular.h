/**
 * Regular grid inherited from ITL_grid.
 * A class which contains information about spatial arrangement
 * of field data points regularly arranged in a Cartesian space.
 * Created on: Dec 9, 2010
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_GRID_REGULAR_H_
#define ITL_GRID_REGULAR_H_

#include "ITL_grid.h"

template <class T>
class ITL_grid_regular: public ITL_grid<T>
{
public:

	/**
	 * Constructor.
	 * @param ndim Dimensionality of the field.
	 */
	ITL_grid_regular( int ndim )
	{
		this->nDim = ndim;
		this->dim = NULL;
	}

	/**
	 * Function for setting the bounds of a grid with no ghost layers.
	 * @param l Pointer to array containing lower grid bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid bounds in continuous space along each dimension.
	 */
	void setBounds( float* l, float* h )
	{
		// set bounding box in physical space
		this->low = new float[this->nDim];
		this->high = new float[this->nDim];
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size to 0 (No pad required)
		this->neighborhoodSize = 0;

		// Set field size and bounding box in terms of grid vertex
		this->dim = new int[this->nDim];
		this->dimWithPad = new int[this->nDim];
		this->lowInt = new int[this->nDim];
		this->highInt = new int[this->nDim];
		this->lowIntWithPad = new int[this->nDim];
		this->highIntWithPad = new int[this->nDim];
		this->lowPad = new int[this->nDim];
		this->highPad = new int[this->nDim];

		// Set the boundary related variables (padding / ghost cell size 0)
		for( int i=0; i<this->nDim; i++ )
		{
			this->lowPad[i] = this->highPad[i] = 0;

			this->lowInt[i] = (int)floor( this->low[i] );
			this->highInt[i] = (int)ceil( this->high[i] );

			this->dim[i] = ( this->highInt[i] - this->lowInt[i] + 1 );
			this->dimWithPad[i] = this->dim[i];

			this->lowIntWithPad[i] = this->lowInt[i];
			this->highIntWithPad[i] = this->highInt[i];
		}

		// Compute number of vertices in the block
		this->nVertices = ITL_util<int>::prod( this->dim, this->nDim );
		this->nVerticesWithPad = this->nVertices;

		#ifdef DEBUG_MODE
		printf( "Field Bounding box:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%f, %f)\t", this->low[i], this->high[i] );
		printf( "\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%d, %d)\t", this->lowInt[i], this->highInt[i] );
		printf( "\n" );
		printf( "Field Dimension:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dim[i] );
		printf( "\n" );
		printf( "Field dimension with pad:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dimWithPad[i] );
		printf( "\n" );
		#endif

	}// end function

	/**
	 * Function for setting the bounds of a grid with ghost layers.
	 * @param l Pointer to array containing lower grid bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layaer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layaer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	void setBounds( float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// set bounding box in physical space
		this->low = new float[this->nDim];
		this->high = new float[this->nDim];
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size and vertices in the neighborhood
		this->neighborhoodSize = neighborhoodsize;

		// Set field size and bounding box in terms of grid vertex
		this->dim = new int[this->nDim];
		this->dimWithPad = new int[this->nDim];
		this->lowInt = new int[this->nDim];
		this->highInt = new int[this->nDim];
		this->lowIntWithPad = new int[this->nDim];
		this->highIntWithPad = new int[this->nDim];
		this->lowPad = new int[this->nDim];
		this->highPad = new int[this->nDim];

		// Set the padding / ghost cell size
		memcpy( this->lowPad, lPad, this->nDim * sizeof(int) );
		memcpy( this->highPad, hPad, this->nDim * sizeof(int) );

		for( int i=0; i<this->nDim; i++ )
		{
			this->lowInt[i] = (int)floor( this->low[i] );
			this->highInt[i] = (int)ceil( this->high[i] );

			this->dim[i] = ( this->highInt[i] - this->lowInt[i] + 1 );
			this->dimWithPad[i] = this->dim[i] + this->highPad[i] + this->lowPad[i];

			this->lowIntWithPad[i] = this->lowInt[i] - this->lowPad[i];
			this->highIntWithPad[i] = this->highInt[i] + this->highPad[i];
		}

		// Compute number of vertices in the block
		this->nVertices = ITL_util<int>::prod( this->dim, this->nDim );
		this->nVerticesWithPad = ITL_util<int>::prod( this->dimWithPad, this->nDim );

		#ifdef DEBUG_MODE
		printf( "Field Bounding box:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%f, %f)\t", this->low[i], this->high[i] );
		printf( "\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%d, %d)\t", this->lowInt[i], this->highInt[i] );
		printf( "\n" );
		printf( "Field Dimension:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dim[i] );
		printf( "\n" );
		printf( "Field dimension with pad:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dimWithPad[i] );
		printf( "\n" );
		#endif

	}// end function

	void setBounds( float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// set bounding box in physical space
		this->low = new float[this->nDim];
		this->high = new float[this->nDim];
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size and vertices in the neighborhood
		this->neighborhoodSizeArray = new int[3];
		this->neighborhoodSizeArray[0] = neighborhoodsizearray[0];
		this->neighborhoodSizeArray[1] = neighborhoodsizearray[1];
		this->neighborhoodSizeArray[2] = neighborhoodsizearray[2];

		// Set field size and bounding box in terms of grid vertex
		this->dim = new int[this->nDim];
		this->dimWithPad = new int[this->nDim];
		this->lowInt = new int[this->nDim];
		this->highInt = new int[this->nDim];
		this->lowIntWithPad = new int[this->nDim];
		this->highIntWithPad = new int[this->nDim];
		this->lowPad = new int[this->nDim];
		this->highPad = new int[this->nDim];

		// Set the padding / ghost cell size
		memcpy( this->lowPad, lPad, this->nDim * sizeof(int) );
		memcpy( this->highPad, hPad, this->nDim * sizeof(int) );

		for( int i=0; i<this->nDim; i++ )
		{
			this->lowInt[i] = (int)floor( this->low[i] );
			this->highInt[i] = (int)ceil( this->high[i] );

			this->dim[i] = ( this->highInt[i] - this->lowInt[i] + 1 );
			this->dimWithPad[i] = this->dim[i] + this->highPad[i] + this->lowPad[i];

			this->lowIntWithPad[i] = this->lowInt[i] - this->lowPad[i];
			this->highIntWithPad[i] = this->highInt[i] + this->highPad[i];
		}

		// Compute number of vertices in the block
		this->nVertices = ITL_util<int>::prod( this->dim, this->nDim );
		this->nVerticesWithPad = ITL_util<int>::prod( this->dimWithPad, this->nDim );

		#ifdef DEBUG_MODE
		printf( "Field Bounding box:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%f, %f)\t", this->low[i], this->high[i] );
		printf( "\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "(%d, %d)\t", this->lowInt[i], this->highInt[i] );
		printf( "\n" );
		printf( "Field Dimension:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dim[i] );
		printf( "\n" );
		printf( "Field dimension with pad:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->dimWithPad[i] );
		printf( "\n" );
		#endif

	}// end function

	/**
	 * Function for 3D spatial index to 1D array index conversion.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 */
	int convert3DIndex( int x, int y, int z )
	{
		return z * this->dimWithPad[0] * this->dimWithPad[1] + y * this->dimWithPad[0] + x;
	}// end function

	/**
	 * Function for 3D+T spatial index to 1D array index conversion.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 * @param t t-coordinate of the spatial point.
	 */
	int convertTimeVarying3DIndex( int x, int y, int z, int t )
	{
		return -1;
	}// end function

	/**
	 * Destructor
	 */
	~ITL_grid_regular()
	{
	}
};

#endif
/* ITL_GRID_REGULAR_H_ */
