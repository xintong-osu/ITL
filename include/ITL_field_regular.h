/**
 * Regular field class inherited from ITL_field.
 * Container for regularly arranged static/time-varying scalar/vector data.
 * Created on: Dec 09, 2010
 * @author Abon
 * @author Teng-Yok
 * @see ITL_field
 */

#ifndef ITL_FIELD_REGULAR_H_
#define ITL_FIELD_REGULAR_H_

#include "ITL_util.h"
#include "ITL_field.h"
#include "ITL_grid_regular.h"
#include "ITL_datastore.h"

template <class T>
class ITL_field_regular: public ITL_field<T>
{
	ITL_grid_regular<T> grid;			/**< Grid associated to the field. */
	ITL_datastore<T> datastore;			/**< Data store associated to the field. */

public:

	/**
	 * Default constructor.
	 */
	ITL_field_regular()
	{
		//this->grid = NULL;
		//this->datastore = NULL;
	}

	ITL_field_regular( const ITL_field_regular<T>& that )
	{
		this->grid = that.grid;
		this->datastore = that.datastore;

	}

	ITL_field_regular( ITL_field_regular<T>& that )
	{
		this->grid = that.grid;
		this->datastore = that.datastore;

	}

	ITL_field_regular& operator= ( const ITL_field_regular<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->grid = that.grid;
			this->datastore = that.datastore;
		}
		// by convention, always return *this
		return *this;
	}

	ITL_field_regular& operator= ( ITL_field_regular<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->grid = that.grid;
			this->datastore = that.datastore;

			/*
			int nDim = that.getNumDim();

			this->grid.init( nDim );

			float l[nDim], h[nDim];
			that.getBounds( l, h );
			this->setBounds( l, h );

			this->datastore.init( that.getDataFull(), that.getSize() );
			*/

		}
		// by convention, always return *this
		return *this;
	}

	/**
	 * Constructor. 
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	ITL_field_regular( int ndim, float* l, float* h )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );

	}// constructor

	/**
	 * Constructor.
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	ITL_field_regular( int ndim, float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsize );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );


	}// Constructor

	/**
	 * Constructor.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsizearray Neighborhood length for each dimension.
	 */
	ITL_field_regular( int ndim, float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );


	}// Constructor

	/**
	 * Constructor.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 */
	ITL_field_regular( T* data, int ndim, int* dim )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds: 0 to N-1
		float* lowEnd = new float[ndim];
		float* highEnd = new float[ndim];
		for( int i=0; i<ndim; i++ )
		{
			lowEnd[i] = 0.0f;
			highEnd[i] = (float)(dim[i]-1);
		}
		this->setBounds( lowEnd, highEnd );

		// Initialize data store
		//this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( data, ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );


	}// constructor 1

	/**
	 * Constructor. 
	 * @param ndim Dimensionality of field.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	ITL_field_regular( T* data, int ndim, float* l, float* h )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		int dim[4];
		getSize( dim );
		datastore.init( data, ITL_util<int>::prod( dim, ndim ) );


	}// constructor

	/**
	 * Constructor.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	ITL_field_regular( T* data, int ndim,
					   float* l, float* h,
					   int* lPad, int* hPad,
					   int neighborhoodsize )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsize );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( data, ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );


	}// Constructor

	/**
	 * Constructor.
	 * @param data pointer to 1D array of elements.
	 * @param ndim Dimensionality of field.
	 * @param dim Length of field along each dimension.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsizearray Neighborhood length for each dimension.
	 */
	ITL_field_regular( T* data, int ndim,
					   float* l, float* h,
					   int* lPad, int* hPad,
					   int* neighborhoodsizearray )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h, lPad, hPad, neighborhoodsizearray );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dimWithPad[0], this->grid->nDim ) );
		int dimWithPad[4];
		getSizeWithPad( dimWithPad );
		datastore.init( data, ITL_util<int>::prod( dimWithPad, grid.getNumDim() ) );


	}// Constructor

	void initialize( int ndim, float *l, float *h )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( ITL_util<int>::prod( &this->grid->dim[0], ndim ) );
		int dim[4];
		getSizeWithPad( dim );
		datastore.init( ITL_util<int>::prod( dim, grid.getNumDim() ) );
	}

	void initialize( T* data, int ndim, float *l, float *h )
	{
		// Initialize grid
		//this->grid = new ITL_grid_regular<T>( ndim );
		grid.init( ndim );

		// set grid bounds
		this->setBounds( l, h );

		// Initialize datastore
		//this->datastore = new ITL_datastore<T>( data, ITL_util<int>::prod( &this->grid->dim[0], ndim ) );
		int dim[4];
		getSize( dim );
		datastore.init( data, ITL_util<int>::prod( dim, grid.getNumDim() ) );

	}



	/**
	 * Pure virtual function for setting the bounds of a field with no ghost layers.
	 * Calls corresponding function of the grid.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 */
	void setBounds( float* l, float* h )
	{	
		// call same method of the grid
		//this->grid->setBounds( l, h );
		grid.setBounds( l, h );

	}// end function

	/**
	 * Pure virtual function for setting the bounds of a field with ghost layers.
	 * Calls corresponding function of the grid.
	 * @param l Pointer to array containing lower grid (associated to the field) bounds in continuous space along each dimension.
	 * @param h Pointer to array containing upper grid (associated to the field) bounds in continuous space along each dimension.
	 * @param lPad Pointer to array containing ghost layer span along each dimension on the lower end.
	 * @param hPad Pointer to array containing ghost layer span along each dimension on the upper end.
	 * @param neighborhoodsize Neighborhood length for each point.
	 */
	void setBounds( float* l, float* h, int* lPad, int* hPad, int neighborhoodsize )
	{
		// call same method of the grid
		//this->grid->setBounds( l, h , lPad, hPad, neighborhoodsize );
		grid.setBounds( l, h , lPad, hPad, neighborhoodsize );

	}// end function

	void setBounds( float* l, float* h, int* lPad, int* hPad, int* neighborhoodsizearray )
	{
		// call same method of the grid
		//this->grid->setBounds( l, h , lPad, hPad, neighborhoodsizearray );
		//float low[4], high[4];
		//int lowpad[4], highpad[4], narray[4];

		//memcpy( low, l, grid.getNumDim() * sizeof(float) );
		//memcpy( high, h, grid.getNumDim() * sizeof(float) );
		//memcpy( narray, neighborhoodsizearray, grid.getNumDim() * sizeof(int) );
		//memcpy( lowpad, lPad, grid.getNumDim() * sizeof(int) );
		//memcpy( highpad, hPad, grid.getNumDim() * sizeof(int) );

		//grid.setBounds( low, high, lowpad, highpad, narray );
		grid.setBounds( l, h, lPad, hPad, neighborhoodsizearray );

	}// end function

	/**
	 * Function for creating a partition or subfield
	 * within the specified bounds.
	 * @param nblocks Number of blocks along each dimension.
	 * @paam subfieldArray A null pointer or a pointer to sub-blocks to be created.  
	 * @param isSharingNeighbor Flag that denotes if two adjacent blocks are sharing a common neighbor. Default value is false.
	 */
	virtual void
	partitionField ( int* nBlock, ITL_field_regular<T>** subfieldArray, bool isSharingNeighbor = true )
	{
		int nDim = this->getNumDim();
		float *blockSize = new float[nDim];
		float* lowSub = new float[nDim];
		float* highSub = new float[nDim];

		// Compute the dimension of each block
		int dim[3];
		this->getSize( dim );
		blockSize[0] = ceil( (dim[0]-1) / (float)nBlock[0] );
		blockSize[1] = ceil( (dim[1]-1) / (float)nBlock[1] );
		blockSize[2] = ceil( (dim[2]-1) / (float)nBlock[2] );
		//blockSize[0] = floor( (this->grid->dim[0]-1) / (float)nBlock[0] );
		//blockSize[1] = floor( (this->grid->dim[1]-1) / (float)nBlock[1] );
		//blockSize[2] = floor( (this->grid->dim[2]-1) / (float)nBlock[2] );


		//#ifdef DEBUG_MODE
		printf( "Dimension of field: %d, %d, %d\n", dim[0], dim[1], dim[2] );
		printf( "Number of blocks (ie, subfields) to be created: %d, %d, %d\n", nBlock[0], nBlock[1], nBlock[2] );
		printf( "Dimension of a subfield: %f, %f, %f\n", blockSize[0], blockSize[1], blockSize[2] );
		//#endif

		// Compute total number of blocks
		int nTotalBlocks = ITL_util<int>::prod( nBlock, nDim );

		// Allocate memory for array of subfields
		assert( (*subfieldArray) != NULL );

		// Partitioning loop
		int blockIndex = 0;
		int lowInt[3], highInt[3];
		this->getBounds( lowInt, highInt );
		for( int z=0; z<nBlock[2]; z++ )
		{
			// Determine [zmin, zmax] for the next block
			lowSub[2] = z * blockSize[2];			
			if( isSharingNeighbor == true )		highSub[2] = std::min( (float)highInt[2], lowSub[2] + blockSize[2] );
						
			for(int y=0; y<nBlock[1]; y++ )
			{
				// Determine [ymin, ymax] for the next block
				lowSub[1] = y * blockSize[1];
				if( isSharingNeighbor == true ) highSub[1] = std::min( (float)highInt[1], lowSub[1] + blockSize[1] );
				
				for( int x=0; x<nBlock[0]; x++)
				{
					// Determine [xmin, xmax] for the next block
					lowSub[0] = x * blockSize[0];
					if( isSharingNeighbor == true )	highSub[0] = std::min( (float)highInt[0], lowSub[0] + blockSize[0] );
					
					// Initialize subfield 
					((*subfieldArray)+blockIndex)->initialize( nDim, lowSub, highSub );//, this->grid->lowPad, this->grid->highPad, this->grid->neighborhoodSize );
					//#ifdef DEBUG_MODE
					printf( "%d: %f %f %f %f %f %f\n", blockIndex, lowSub[0], lowSub[1], lowSub[2], highSub[0], highSub[1], highSub[2] );
					//#endif

					// load (copy from datastore of full field) Data to the subfield
					cout << "a" << endl;
					int nSubPoint = ( highSub[0] - lowSub[0] + 1 ) * ( highSub[1] - lowSub[1] + 1 ) * ( highSub[2] - lowSub[2] + 1 );
					T* subData = new T[nSubPoint];
					this->getDataBetween( lowSub, highSub, subData );
					((*subfieldArray)+blockIndex)->setDataFull( subData, nSubPoint );
					delete [] subData;
					cout << "b" << endl;

					//#ifdef DEBUG_MODE
					float m = ITL_util<SCALAR>::Min( (SCALAR*)((*subfieldArray)+blockIndex)->getDataFull(), ((*subfieldArray)+blockIndex)->getSize() );
					float M = ITL_util<SCALAR>::Max( (SCALAR*)((*subfieldArray)+blockIndex)->getDataFull(), ((*subfieldArray)+blockIndex)->getSize() );
					cout << m << " " << M << endl;
					//#endif

					// Increment to the next block
					blockIndex ++;
				}
			}
		}

		delete [] lowSub;
		delete [] highSub;
	}// end function


	/**
	 * Function for creating a partition or subfield
	 * within the specified bounds.
	 * @param low Lower bound along each dimension.
	 * @param high Higher bound along each dimension.
	 */
	/*
	ITL_field<T>* createSubField( float* low, float* high )
	{
		// Allocate memory for 3D block
		ITL_field_regular<T>* newField = new ITL_field_regular<T>( this->getDataBetween( low, high ), this->grid->nDim, this->grid->dim );

		return newField;

	}// end function
	*/

	/**
	 * Data accessor function type 1.
	 * @param id 1D index of data array
	 * @return field value.
	 */
	virtual T getDataAt( int id )
	{
		// Get data
		//return this->datastore->array[id];
		return datastore.getDataAt(id);
	}// end function

	/**
	 * Data accessor function type 2.
	 * @param x x-coordinate of the spatial point.
	 * @param y y-coordinate of the spatial point.
	 * @param z z-coordinate of the spatial point.
	 * @return field value.
	 */
	virtual T getDataAt( int x, int y, int z )
	{
		// Convert index to 1D
		//int index1d = this->grid->convert3DIndex( x, y, z );
		int index1d = grid.convert3DIndex( x, y, z );

		// Get data
		//return this->datastore->array[index1d];
		return datastore.getDataAt( index1d );
	}// end function

	/**
	 * Data accessor function type 3.
	 * Returns chunk of data within specified bound.
	 * @param lowBoundary Lower bound along each dimension.
	 * @param highBoundary Higher bound along each dimension.
	 * @return pointer to field value.
	 */
	virtual void
	getDataBetween( float* lowBoundary, float* highBoundary, T* retData )
	{
		// ADD-BY-LEETEN 02/13/2012-BEGIN
		#ifdef	WIN32
			int* lowBoundaryInt = new int[this->grid->nDim];
			int* highBoundaryInt = new int[this->grid->nDim];
		#else	// #ifdef	WIN32
		// ADD-BY-LEETEN 02/13/2012-END

			//int lowBoundaryInt[this->grid->nDim];
			//int highBoundaryInt[this->grid->nDim];
			int lowBoundaryInt[grid.getNumDim()];
			int highBoundaryInt[grid.getNumDim()];

		#endif	// #ifdef	WIN32		// ADD-BY-LEETEN 02/13/2012

		//for( int i=0; i<this->grid->nDim; i++ )
		for( int i=0; i<grid.getNumDim(); i++ )
		{
			lowBoundaryInt[i] = (int)floor( lowBoundary[i] );
			highBoundaryInt[i] = (int)ceil( highBoundary[i] );
		}

		this->getDataBetween( lowBoundaryInt, highBoundaryInt, retData );
		
		// ADD-BY-LEETEN 02/13/2012-BEGIN
		#ifdef	WIN32
		delete [] lowBoundaryInt;
		delete [] highBoundaryInt;
		#endif	// #ifdef	WIN32
		// ADD-BY-LEETEN 02/13/2012-END
	}// end function

	/**
	 * Data accessor function type 4.
	 * Returns chunk of data within specified bound.
	 * @param lowBoundary Lower bound along each dimension.
	 * @param highBoundary Higher bound along each dimension.
	 * @return pointer to field value.
	 */
	virtual void
	getDataBetween( int* lowBoundary, int* highBoundary, T* retData )
	{
		int dimLength[4];
		int blockLowInt[3], blockHighInt[3], blockSize[3];
		int blockOffset = 0, localOffset = 0;

		// Count number of vertices requested
		ITL_util<int>::subtractArrays( highBoundary, lowBoundary, dimLength, grid.getNumDim() );
		ITL_util<int>::addArrayScalar( dimLength, 1, grid.getNumDim() );
		int nV = ITL_util<int>::prod( dimLength, grid.getNumDim() );

		// Get bounds of the grid
		grid.getBounds( blockLowInt, blockHighInt );

		#ifdef DEBUG_MODE
		grid.getSize( blockSize );
		fprintf( stderr, "block boundary: %d %d %d, %d %d %d\n", blockLowInt[0], blockLowInt[1], blockLowInt[2],
																 blockHighInt[0], blockHighInt[1], blockHighInt[2] );
		fprintf( stderr, "%d %d %d\n", blockSize[0], blockSize[1], blockSize[2] );
		#endif

		for( int z=0; z<dimLength[2]; z++ )
		{
			for( int y=0; y<dimLength[1]; y++ )
			{
				// DELETED-BY-ABON START: 08/29/12
				//blockOffset = grid.convert3DIndex( lowBoundary[0],
				//								   lowBoundary[1] + y,
				//								   lowBoundary[2] + z );
				// DELETED-BY-ABON END: 08/29/12

				// LowBoundary and highBoundary are in global space.
				// The need to be converted to the block's local space by
				// translating relative to the current block's start point.
				blockOffset = grid.convert3DIndex( lowBoundary[0]  - blockLowInt[0],
												   lowBoundary[1]  - blockLowInt[1] + y,
												   lowBoundary[2]  - blockLowInt[2] + z );
				//blockOffset = ( lowBoundary[2] - blockLowInt[2] + z ) * blockSize[0] * blockSize[1]
				 //             + ( lowBoundary[1] - blockLowInt[1] + y ) * blockSize[0]
				  //            + ( lowBoundary[0]  - blockLowInt[0] );

				//fprintf( stderr, "%d %d %d %d\n", y, z, blockOffset, localOffset );
				//if( blockOffset + dimLength[0] > blockSize[0]*blockSize[1]*blockSize[2] )
				//{

				//	fprintf( stderr, "%d %d %d\n", lowBoundary[0]  - blockLowInt[0],
				//								   lowBoundary[1] - blockLowInt[1] + y,
				//								   lowBoundary[2] - blockLowInt[2] + z );
				//	fprintf( stderr, "%d %d %d\n", dimLength[0], dimLength[2], dimLength[2] );
				//	fprintf( stderr, "%d %d %d\n", blockSize[0], blockSize[2], blockSize[2] );
				//	fprintf( stderr, "%d %d %d %d\n", y, z, blockOffset, localOffset );
				//	fprintf( stderr, "%d %d\n", blockOffset + dimLength[0], blockSize[0]*blockSize[1]*blockSize[2]);
				//}
				//fprintf( stderr, "%d %d %d %d\n", y, z, blockOffset, localOffset );

				memcpy( (T*)(retData + localOffset ) , (T*)(getDataFull()+blockOffset ) , dimLength[0] * sizeof( T ) );
				localOffset += dimLength[0];
			}
		}


	}// end function
	
	/**
	 * Data accessor function type 5.
	 * Returns entire field data.
	 * @return pointer to field value.
	 */
	virtual T* getDataFull()
	{
		//return this->datastore->array;
		return datastore.getData();
	}// end function

	virtual int
	getSize()
	{
		return grid.getSize();
	}

	virtual void
	getSize( int* d )
	{
		int dim[4];
		grid.getSize( dim );

		memcpy( d, dim, sizeof(int)*grid.getNumDim() );
	}

	virtual void
	getSizeWithPad( int* d )
	{
		int dimWPad[4];
		grid.getSizeWithPad( dimWPad );

		memcpy( d, dimWPad, sizeof(int)*grid.getNumDim() );
	}

	virtual int
	getNumDim()
	{
		return grid.getNumDim();
	}

	virtual void
	getBounds( float* l, float* h )
	{
		float low[4];
		float high[4];
		grid.getBounds( low, high );

		memcpy( l, low, sizeof(float)*grid.getNumDim() );
		memcpy( h, high, sizeof(float)*grid.getNumDim() );
	}

	virtual void
	getBounds( int* l, int* h )
	{
		int low[4];
		int high[4];
		grid.getBounds( low, high );

		memcpy( l, low, sizeof(int)*grid.getNumDim() );
		memcpy( h, high, sizeof(int)*grid.getNumDim() );
	}

	/**
	 * Field bounds accessor function.
	 */
	virtual void
	getBounds( int *limits )
	{
		int lim[6];
		grid.getBounds( lim );

		limits[0] = lim[0];
		limits[1] = lim[1];
		limits[2] = lim[2];
		limits[3] = lim[3];
		limits[4] = lim[4];
		limits[5] = lim[5];

	}

	virtual void
	getPadSize( int* lowpad, int* highpad )
	{
		int lowPad[4];
		int highPad[4];

		grid.getPadSize( lowPad, highPad );

		memcpy( lowpad, lowPad, sizeof(int)*grid.getNumDim() );
		memcpy( highpad, highPad, sizeof(int)*grid.getNumDim() );

	}

	virtual int
	getNeighborhoodSize()
	{
		return grid.getNeighborhoodSize();
	}

	virtual void
	getNeighborhoodSize( int* neighborsize )
	{
		int neighborSize[4];

		grid.getNeighborhoodSize( neighborSize );

		memcpy( neighborsize, neighborSize, sizeof(int)*grid.getNumDim() );

	}



	/**
	 * Data mutator function type 1.
	 * Sets data at a particular array location.
	 * @param id 1D index of data array
	 * @param data field value.
	 */
	virtual void setDataAt( int id, T data )
	{
		// Set data
		//this->datastore->array[id] = data;
		datastore.setDataAt( data, id );
	}// end function

	/**
	 * Data mutator function type 2.
	 */
	virtual void setDataAt( int x, int y, int z, T data )
	{
		// Convert index to 1D
		//int index1d = this->grid->convert3DIndex( x, y, z );
		int index1d = grid.convert3DIndex( x, y, z );

		// Set data
		//this->datastore->array[index1d] = data;
		datastore.setDataAt( data, index1d );

	}// end function

	/**
	 * Data mutator function type 3.
	 */
	virtual void setDataBetween( T* data )
	{
		// Not implemented yet
	}// end function

	/**
	 * Data mutator function type 4.
	 * Sets data for entire field (shallow copy).
	 * @param data pointer to data
	 */
	virtual void
	setDataFull( T* data )
	{
		//this->datastore->array = data;
		datastore.setDataFull( data );
	}// end function

	/**
	 * Data mutator function type 4.
	 * Sets data for entire field (deep copy).
	 * @param data pointer to data
	 */
	virtual void
	setDataFull( T* data, int ndata )
	{
		//this->datastore->array = data;
		datastore.setDataFull( data, ndata );
	}// end function


	int
	convert3DIndex( int x, int y, int z )
	{
		return grid.convert3DIndex( x, y, z );
	}

	/**
	 * Destructor.
	 */
	virtual
	~ITL_field_regular()
	{
		//delete grid;
		//delete datastore;
	}
};

#endif
/* ITL_FIELD_REGULAR_H_ */
