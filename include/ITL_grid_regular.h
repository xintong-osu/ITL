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
	int nDim;			/**< Number of dimensions */
	int dim[4];			/**< Array for length of each dimension */
	int dimWithPad[4];		/**< Array for length of each dimension along with ghost layers (if any) */
	int nVertices;			/**< Number of vertices in Cartesian space. */
	int nVerticesWithPad;		/**< Number of vertices along with ghost layers (if any) in Cartesian space. */
	float low[4];			/**< Lower bound of field in continuous Cartesian space. Limit may not coincide with a grid vertex */
	float high[4];			/**< Upper bound of field in continuous Cartesian space. Limit may not coincide with a grid vertex */
	int lowInt[4];			/**< Lower bound of field in discrete Cartesian space (nearest grid vertex containing the lower bound). */
	int highInt[4];			/**< Upper bound of field in discrete Cartesian space (nearest grid vertex containing the upper bound). */
	int lowIntWithPad[4];		/**< Lower bound of field in discrete Cartesian space along with ghost layers (if any) */
	int highIntWithPad[4];		/**< Upper bound of field in discrete Cartesian space along with ghost layers (if any) */

	int neighborhoodSize;   	/**< Neighborhood length for each point. Helps to determine span of ghost layers */
	int lowPad[4];			/**< Span of ghost layers beyond lower end. */
	int highPad[4];			/**< Span of ghost layers beyond upper end. */
	int neighborhoodSizeArray[4];

public:

	ITL_grid_regular()
	{
	}

	/**
	 * Constructor.
	 * @param ndim Dimensionality of the field.
	 */
	ITL_grid_regular( int ndim )
	{
		nDim = ndim;
		//this->dim = NULL;
		//this->low = NULL;
		//this->high = NULL;
		//this->lowInt = NULL;
		//this->highInt = NULL;
		//this->lowPad = NULL;
		//this->highPad = NULL;
		//this->lowIntWithPad = NULL;
		//this->highIntWithPad = NULL;
	}

	ITL_grid_regular( const ITL_grid_regular<T>& that )
	{
		this->nDim = that.nDim;
		this->nVertices = that.nVertices;
		this->nVerticesWithPad = that.nVerticesWithPad;
		this->neighborhoodSize = that.neighborhoodSize;

		for( int i=0; i<nDim; i++ )
		{
			this->dim[i] = that.dim[i];
			this->dimWithPad[i] = that.dimWithPad[i];

			this->low[i] = that.low[i];
			this->high[i] = that.high[i];
			this->lowInt[i] = that.lowInt[i];
			this->highInt[i] = that.highInt[i];
			this->lowPad[i] = that.lowPad[i];
			this->highPad[i] = that.highPad[i];
			this->lowIntWithPad[i] = that.lowIntWithPad[i];
			this->highIntWithPad[i] = that.highIntWithPad[i];

			this->neighborhoodSizeArray[i] = that.neighborhoodSizeArray[i];
			this->highIntWithPad[i] = that.highIntWithPad[i];

		}

		//float l[nDim], h[nDim];
		//int lPad[nDim], hPad[nDim];
		//int neighborhoodsizearray[nDim];

		//that.getBounds( l, h );
		//that.getPadSize( lPad, hPad );
		//that.getNeighborSize( neighborhoodsizearray );

		//this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

	}

	ITL_grid_regular( ITL_grid_regular<T>& that )
	{
		this->nDim = that.nDim;
		this->nVertices = that.nVertices;
		this->nVerticesWithPad = that.nVerticesWithPad;
		this->neighborhoodSize = that.neighborhoodSize;

		for( int i=0; i<nDim; i++ )
		{
			this->dim[i] = that.dim[i];
			this->dimWithPad[i] = that.dimWithPad[i];

			this->low[i] = that.low[i];
			this->high[i] = that.high[i];
			this->lowInt[i] = that.lowInt[i];
			this->highInt[i] = that.highInt[i];
			this->lowPad[i] = that.lowPad[i];
			this->highPad[i] = that.highPad[i];
			this->lowIntWithPad[i] = that.lowIntWithPad[i];
			this->highIntWithPad[i] = that.highIntWithPad[i];

			this->neighborhoodSizeArray[i] = that.neighborhoodSizeArray[i];
			this->highIntWithPad[i] = that.highIntWithPad[i];

		}

		/*
		this->nDim = that.getNumDim();

		float l[nDim], h[nDim];
		int lPad[nDim], hPad[nDim];
		int neighborhoodsizearray[nDim];

		that.getBounds( l, h );
		that.getPadSize( lPad, hPad );
		that.getNeighborSize( neighborhoodsizearray );

		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );
		*/

	}

	ITL_grid_regular& operator= ( const ITL_grid_regular<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->nDim = that.nDim;
			this->nVertices = that.nVertices;
			this->nVerticesWithPad = that.nVerticesWithPad;
			this->neighborhoodSize = that.neighborhoodSize;

			for( int i=0; i<nDim; i++ )
			{
				this->dim[i] = that.dim[i];
				this->dimWithPad[i] = that.dimWithPad[i];

				this->low[i] = that.low[i];
				this->high[i] = that.high[i];
				this->lowInt[i] = that.lowInt[i];
				this->highInt[i] = that.highInt[i];
				this->lowPad[i] = that.lowPad[i];
				this->highPad[i] = that.highPad[i];
				this->lowIntWithPad[i] = that.lowIntWithPad[i];
				this->highIntWithPad[i] = that.highIntWithPad[i];

				this->neighborhoodSizeArray[i] = that.neighborhoodSizeArray[i];
				this->highIntWithPad[i] = that.highIntWithPad[i];

			}


			/*
			this->nDim = that.getNumDim();

			float l[nDim], h[nDim];
			int lPad[nDim], hPad[nDim];
			int neighborhoodsizearray[nDim];

			that.getBounds( l, h );
			that.getPadSize( lPad, hPad );
			that.getNeighborSize( neighborhoodsizearray );

			this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );
			*/

		}
		// by convention, always return *this
		return *this;
	}

	ITL_grid_regular& operator= ( ITL_grid_regular<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->nDim = that.nDim;
			this->nVertices = that.nVertices;
			this->nVerticesWithPad = that.nVerticesWithPad;
			this->neighborhoodSize = that.neighborhoodSize;

			for( int i=0; i<nDim; i++ )
			{
				this->dim[i] = that.dim[i];
				this->dimWithPad[i] = that.dimWithPad[i];

				this->low[i] = that.low[i];
				this->high[i] = that.high[i];
				this->lowInt[i] = that.lowInt[i];
				this->highInt[i] = that.highInt[i];
				this->lowPad[i] = that.lowPad[i];
				this->highPad[i] = that.highPad[i];
				this->lowIntWithPad[i] = that.lowIntWithPad[i];
				this->highIntWithPad[i] = that.highIntWithPad[i];

				this->neighborhoodSizeArray[i] = that.neighborhoodSizeArray[i];
				this->highIntWithPad[i] = that.highIntWithPad[i];

			}

			/*
			this->nDim = that.getNumDim();

			float l[nDim], h[nDim];
			int lPad[nDim], hPad[nDim];
			int neighborhoodsizearray[nDim];

			that.getBounds( l, h );
			that.getPadSize( lPad, hPad );
			that.getNeighborSize( neighborhoodsizearray );

			this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );
			*/

		}
		// by convention, always return *this
		return *this;
	}

	void
	init( int ndim )
	{
		nDim = ndim;
	}

	/**
	 * Function for setting the bounds of a grid with no ghost layers.
	 * @param l Pointer to array containing lower grid bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid bounds in continuous space along each dimension.
	 */
	void setBounds( float* l, float* h )
	{
		// set bounding box in physical space
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size to 0 (No pad required)
		this->neighborhoodSize = 0;
		this->neighborhoodSizeArray[0] = this->neighborhoodSizeArray[1] = this->neighborhoodSizeArray[2] = 0;

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
		//this->low = new float[this->nDim];
		//this->high = new float[this->nDim];
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size and vertices in the neighborhood
		this->neighborhoodSize = neighborhoodsize;
		//this->neighborhoodSizeArray = new int[3];
		this->neighborhoodSizeArray[0] = this->neighborhoodSizeArray[1] = this->neighborhoodSizeArray[2] = this->neighborhoodSize;   

		// Set field size and bounding box in terms of grid vertex
		//this->dim = new int[this->nDim];
		//this->dimWithPad = new int[this->nDim];
		//this->lowInt = new int[this->nDim];
		//this->highInt = new int[this->nDim];
		//this->lowIntWithPad = new int[this->nDim];
		//this->highIntWithPad = new int[this->nDim];
		//this->lowPad = new int[this->nDim];
		//this->highPad = new int[this->nDim];

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
		memcpy( this->low, l, this->nDim * sizeof(float) );
		memcpy( this->high, h, this->nDim * sizeof(float) );

		// Set the neighborhood size and vertices in the neighborhood
		memcpy( this->neighborhoodSizeArray, neighborhoodsizearray, this->nDim * sizeof(int) );

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
			printf( "(%g, %g)\t", this->low[i], this->high[i] );
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
		printf( "Neighborhood size:\n" );
		for( int i=0; i<this->nDim; i++ )
			printf( "%d\t", this->neighborhoodSizeArray[i] );
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


	virtual int
	getSize()
	{
		return nVertices;
	}

	virtual void
	getSize( int* d )
	{
		memcpy( d, dim, sizeof(int)*nDim );
	}

	virtual void
	getSizeWithPad( int* d )
	{
		memcpy( d, dimWithPad, sizeof(int)*nDim );
	}

	virtual void
	getBounds( float* l, float* h )
	{
		memcpy( l, low, sizeof(float)*nDim );
		memcpy( h, high, sizeof(float)*nDim );
	}

	virtual void
	getBounds( int* l, int* h )
	{
		memcpy( l, lowInt, sizeof(int)*nDim );
		memcpy( h, highInt, sizeof(int)*nDim );
	}

	virtual void
	getBounds( int *limits )
	{
		limits[0] = lowInt[0];
		limits[1] = lowInt[1];
		limits[2] = lowInt[2];
		limits[3] = highInt[0];
		limits[4] = highInt[1];
		limits[5] = highInt[2];

	}

	virtual void
	getPadSize( int* lowpad, int* highpad )
	{
		memcpy( lowpad, lowPad, sizeof(int)*getNumDim() );
		memcpy( highpad, highPad, sizeof(int)*getNumDim() );

	}

	virtual int
	getNeighborhoodSize()
	{
		return neighborhoodSize;
	}


	virtual void
	getNeighborhoodSize( int* neighborsize )
	{
		memcpy( neighborsize, neighborhoodSizeArray, sizeof(int)*getNumDim() );
	}


	virtual int
	getNumDim()
	{
		return nDim;
	}



	/**
	 * Destructor
	 */
	virtual
	~ITL_grid_regular()
	{
		//if( this->dim != NULL )		delete this->dim;
		//if( this->dimWithPad != NULL ) 	delete this->dimWithPad;
		//if( this->low != NULL )		delete this->low;
		//if( this->high != NULL )	delete this->high;
		//if( this->lowInt != NULL )	delete this->lowInt;
		//if( this->highInt != NULL )	delete this->highInt;
		//if( this->lowIntWithPad != NULL )	delete  this->lowIntWithPad;
		//if( this->highIntWithPad != NULL )	delete  this->highIntWithPad;
		//if( this->lowPad != NULL )	delete  this->lowPad;
		//if( this->highPad != NULL )	delete  this->highPad;
	}// end destructor
};

#endif
/* ITL_GRID_REGULAR_H_ */
