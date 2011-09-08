/**
 * Template-based Utility class for File I/O.
 * Implements static functions for serial and parallel file I/O.
 * Created on: Nov 23, 2010
 * @author Abon
 * @author Teng-Yok
 */

#ifndef ITL_IOUTIL_H_
#define ITL_IOUTIL_H_

#include <mpi.h>
#include "ITL_header.h"
#include "ITL_util.h"

template <class T>
class ITL_ioutil
{
public:

	/**
	 * Serial file reader.
	 * Reads field stored as binary file in vec format.
	 * @param fileName File name.
	 * @param nDim Dimensionality of the field.
	 * @param dim Length of field along each dimension (can be a zero filled array)
	 * @return pointer to data
	 */
	static T* readFieldBinarySerial( const char* fileName, int nDim, int* dim )
	{
		int nel = 0;

		// Open file
		FILE* dataFile = fopen( fileName, "rb" );
		assert( dataFile != NULL );

		// Read header
		for( int i=0; i<nDim; i++ )
			fread( &dim[i], sizeof(int), 1, dataFile );
				
		// Compute total number of elements
		nel = ITL_util<int>::prod( dim, nDim );

		// Allocate memory for data to be read
		T* array = new T[nel];

		// Read data to 1D array
		fread( &array[0], sizeof(T), nel, dataFile );

		#ifdef DEBUG_MODE
		if( nDim == 3 )
			printf( "NX: %d, NY: %d, NZ: %d\nTotal number of elements: %d\n", dim[0], dim[1], dim[2], nel );
		printf( "%d values read from file\n", nel );
		#endif

		// close file
		fclose( dataFile );

		// return data array
		return array;

	}// end function

	/**
	 * Serial file writer.
	 * Stores field binary file in vec format.
	 * @param data pointer to data
	 * @param fileName File name.
	 * @param dim Length of field along each dimension (can be a zero filled array)
	 * @param nFieldDim Dimensionality of the field
	 */
	static void writeFieldBinarySerial( T* data, const char* fileName, int* dim, int nFieldDim )
	{
		// Compute total number of scalar values to write
		int nel = ITL_util<int>::prod( dim, nFieldDim );

		// Open file for saving entropy field
		FILE* outFile = fopen( fileName, "wb" );
		assert( fileName != NULL );

		// Write scalar field dimensions to file header
		for( int i=0; i<nFieldDim; i++ )
			fwrite( &dim[i], sizeof(int), 1, outFile );

		// Write data
		assert( data != NULL );
		fwrite( &data[0], sizeof(T), nel, outFile );

		// close file
		fclose( outFile );
	}

	/**
	 * Parallel file reader.
	 * Reads field stored as binary file in vec format.
	 * @param fileName File name.
	 * @param nDim Dimensionality of the field.
	 * @param dataDim Length of field along each dimension (can be a zero filled array).
	 * @param blockDim Length of the block along each dimension to be read by this instance
	 * @param nBlocks Number ob blocks along each dimension into which the entire field has
	 * been subdivided (can be a zero filled array).
	 * of the program (can be a zero filled array).
	 * @param low Pointer to array containing lower grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param high Pointer to array containing upper grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param lowPad Pointer to array containing ghost layer span along each dimension on the
	 * lower end (can be a zero filled array).
	 * @param highPad Pointer to array containing ghost layer span along each dimension on the
	 * upper end (can be a zero filled array).
	 * @param neighborhoodSize Neighborhood length for each point.
	 * @param myId Rank of this processor
	 * @param nProcs Total number of processes
	 * @return pointer to portion of data
	 */
	static T* readFieldBinaryParallel( char* fileName, int nDim, int* dataDim,
									   int* blockDim, int* nBlocks, int* blockId,
									   int* low, int* high,
									   int* lowPad, int* highPad, int neighborhoodSize,
									   int myId, int nProcs )
	{
		int nel = 0, nelBlock = 0, nelBlockWithPad = 0;
		int file_open_error;
		MPI_Status status;
		MPI_Offset offset;

		int* paddedLow = new int[nDim];
		int* paddedHigh = new int[nDim];
		int* paddedBlockDim = new int[nDim];

		// Open file
		MPI_File dataFile;
		file_open_error = MPI_File_open( MPI_COMM_WORLD, fileName,
				                  	  	 MPI_MODE_RDONLY, MPI_INFO_NULL, &dataFile );
		assert( dataFile != NULL && file_open_error == MPI_SUCCESS );

		// Set view ** Not necessary

		// read header to know the dimension of data
		char buffer[4];
		for( int i=0; i<nDim; i++ )
		{
			MPI_File_read_at( dataFile, i*4, buffer, 4, MPI_CHAR, &status );
			dataDim[i] = *((int*)buffer);
		}

		// Convert 1D myId to 3D block ID
		int quotient = myId;
		for( int i=0; i<nDim; i++ )
		{
			blockId[i] = quotient % nBlocks[i];
			quotient =  quotient / nBlocks[i];
		}

		// Compute size and bounding box of each block
		assert( ITL_util<int>::prod( nBlocks, nDim ) == nProcs );
		for( int i=0; i<nDim; i++ )
		{
			blockDim[i] = (int) floor( (float)dataDim[i] / nBlocks[i] );
			low[i] = blockId[i] * blockDim[i];

			// Make necessary adjustments if the last processor gets a block
			// of different size
			if( dataDim[i] % nBlocks[i] != 0 && blockId[i] == nBlocks[i]-1 )
				blockDim[i] = blockDim[i] + dataDim[i] % nBlocks[i];

			high[i] = ( low[i] + blockDim[i]-1 );
			paddedLow[i] = ITL_util<int>::clamp( low[i] - neighborhoodSize, 0, dataDim[i]-1 );
			paddedHigh[i] = ITL_util<int>::clamp( high[i] + neighborhoodSize, 0, dataDim[i]-1 );
			paddedBlockDim[i] = paddedHigh[i] - paddedLow[i] + 1;

			lowPad[i] = low[i] - paddedLow[i];
			highPad[i] = paddedHigh[i] - high[i];

		}// end for

		// Compute number of elements in data and in this block
		nel = ITL_util<int>::prod( dataDim, nDim );
		nelBlock = ITL_util<int>::prod( blockDim, nDim );
		nelBlockWithPad = ITL_util<int>::prod( paddedBlockDim, nDim );

		#ifdef DEBUG_MODE
		if( nDim == 3 )
		{
			printf( "%d: NX: %d, NY: %d, NZ: %d\nTotal number of elements in data: %d\n", myId, dataDim[0], dataDim[1], dataDim[2], nel );
			printf( "%d: Block ID: %d %d %d\n", myId, blockId[0], blockId[1], blockId[2] );
			printf( "%d: Block Start point: %d %d %d\n", myId, low[0], low[1], low[2] );
			printf( "%d: Block End point: %d %d %d\n", myId, high[0], high[1], high[2] );
			printf( "%d: Block Start point (with pad): %d %d %d\n", myId, paddedLow[0], paddedLow[1], paddedLow[2] );
			printf( "%d: Block End point (with pad): %d %d %d\n", myId, paddedHigh[0], paddedHigh[1], paddedHigh[2] );
			printf( "%d: Block Size (with pad): %d %d %d\n", myId, paddedBlockDim[0], paddedBlockDim[1], paddedBlockDim[2] );
			printf( "%d: Lower side Pad Size: %d %d %d\n", myId, lowPad[0], lowPad[1], lowPad[2] );
			printf( "%d: Higher side Pad Size: %d %d %d\n", myId, highPad[0], highPad[1], highPad[2] );
		}
		#endif

		// Allocate memory for data to be read
		T* array = new T[nelBlockWithPad];

		// Read data to 1D array
		int arrayIndex = 0;
		int fileIndex = 0;

		for( int z=paddedLow[2]; z<=paddedHigh[2]; z++ )
		{
			for( int y=paddedLow[1]; y<=paddedHigh[1]; y++ )
			{
				fileIndex = z * dataDim[0] * dataDim[1] + y * dataDim[0] + paddedLow[0];
				offset = nDim * sizeof(int) + fileIndex * sizeof(T);

				MPI_File_read_at( dataFile, offset, array+arrayIndex, sizeof(T)*paddedBlockDim[0], MPI_CHAR, &status );

				// Increment array Index
				arrayIndex += paddedBlockDim[0];

				/*
				// The following section reads one element at a time
				// Commented out by ABON: 12-12-10
				for( int x=paddedLow[0]; x<=paddedHigh[0]; x++ )
				{
					fileIndex = z * dataDim[0] * dataDim[1] + y * dataDim[0] + x;
					offset = nDim * sizeof(int) + fileIndex * sizeof(T);

					// Read Data
					MPI_File_read_at( dataFile, offset, array+arrayIndex, sizeof(T), MPI_CHAR, &status );

					// Increment array Index
					arrayIndex += 1;
				}
				*/

			}
		}

		#ifdef DEBUG_MODE
		printf( "%d: %d values read from file\n", myId, nelBlockWithPad );
		#endif

		// close file
		MPI_File_close( &dataFile );

		// Free Resources
		delete [] paddedLow;
		delete [] paddedHigh;
		delete [] paddedBlockDim;

		// return data array
		return array;

	}// end function

	/**
	 * Parallel file writer.
	 * Stores field as binary file in vec format.
	 * @param data pointer to portion of data
	 * @param fileName File name.
	 * @param dataDim Length of field along each dimension (can be a zero filled array).
	 * @param blockDim Length of the block along each dimension to be read by this instance
	 * of the program (can be a zero filled array).
	 * @param low Pointer to array containing lower grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param nFieldDim Dimensionality of the field.
	 * @param myId Rank of this processor
	 * @param nProcs Total number of processes
	 */
	static void writeFieldBinaryParallel( T* data, char* fileName,
			                              int* dataDim, int* blockDim,
			                              float* low, int nFieldDim,
			                              int myId, int nProcs )
	{
		int file_open_error, file_allocate_error, file_write_error;
		MPI_Status status;
		MPI_Offset offset;

		// Compute total number of values to write
		int nel = ITL_util<int>::prod( dataDim, nFieldDim );

		// Open file
		MPI_File dataFile;
		file_open_error = MPI_File_open( MPI_COMM_WORLD, fileName,
				                  	  	 MPI_MODE_WRONLY | MPI_MODE_CREATE,
				                  	  	 MPI_INFO_NULL, &dataFile );
		assert( dataFile != NULL );
		assert( file_open_error == MPI_SUCCESS );

		// Pre-allocate space for file ( All processes need to allocate same amount of space )
		file_allocate_error = MPI_File_preallocate( dataFile, nFieldDim*sizeof(int) + nel * sizeof(T) );
		assert( file_allocate_error == MPI_SUCCESS );

		// Write scalar field dimensions to file header ( All processes )
		for( int i=0; i<nFieldDim; i++ )
		{
			file_write_error = MPI_File_write_at_all( dataFile, i*sizeof(int), dataDim+i, sizeof(int), MPI_CHAR, &status );
			assert( file_write_error == MPI_SUCCESS );
		}

		// Write data of this block at appropriate position
		assert( data != NULL );
		offset = nFieldDim*sizeof(int);
		int arrayIndex = 0;
		int fileIndex = 0;
		int* lowInt = new int[nFieldDim];
		for( int i=0; i<nFieldDim; i++ )
			lowInt[i] = (int)floor( low[i] );

		for( int z=lowInt[2]; z<(lowInt[2]+blockDim[2]); z++ )
		{
			for( int y=lowInt[1]; y<(lowInt[1]+blockDim[1]); y++ )
			{
				for( int x=lowInt[0]; x<(lowInt[0]+blockDim[0]); x++ )
				{
					fileIndex = z * dataDim[0] * dataDim[1] + y * dataDim[0] + x;
					offset = nFieldDim * sizeof(int) + fileIndex * sizeof(T);

					file_write_error = MPI_File_write_at( dataFile, offset, data+arrayIndex, sizeof(T), MPI_CHAR, &status );
					assert( file_write_error == MPI_SUCCESS );

					// Increment array Index
					arrayIndex += 1;
				}
			}
		}

		// close file
		MPI_File_close( &dataFile );

	}// end function

	/**
	 * Modified version of parallel file reader.
	 * Stores field as binary file in vec format. Uses advanced MPI-2 features like
	 * subarray, file view etc.
	 * @param fileName File name.
	 * @param nDim Number of dimensions of data.
	 * @param dataDim Length of field along each dimension (can be a zero filled array).
	 * @param blockDim Length of the block along each dimension to be read by this instance
	 * of the program (can be a zero filled array).
	 * @param nBlocks Pointer to array containing number of blocks along each dimension.
	 * @param blockId Pointer to array containing serial number of block along each dimension.
	 * @param low Pointer to array containing lower grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param high Pointer to array containing higher grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param lowPad Pointer to array containing ghost layer span along each dimension on the
	 * lower end (can be a zero filled array).
	 * @param highPad Pointer to array containing ghost layer span along each dimension on the
	 * upper end (can be a zero filled array).
	 * @param neighborhoodSize Neighborhood length for each point.
	 * @param myId Rank of this processor
	 * @param nProcs Total number of processes
	 * @return pointer to portion of data
	 */
	static T* readFieldBinaryParallel2( const char* fileName, int nDim, int* dataDim,
								int* blockDim,  int* nBlocks, int* blockId,
								int* low, int* high,
								int* lowPad, int* highPad, int neighborhoodSize,
								int myId, int nProcs )
	{
		MPI_File dataFile;
		MPI_Datatype eType, fileType;
		int fileOpenError, fileReadError, subarrayError, viewSetError, commitError;
		MPI_Status status;
		MPI_Offset offset;
		int* paddedLow = new int[nDim];
		int* paddedHigh = new int[nDim];
		int* paddedBlockDim = new int[nDim];

		// Create subarray for file header
		int *headerSize = new int[1];  headerSize[0] = 3;
		int *headerStart = new int[1]; headerStart[0] = 0;
		subarrayError = MPI_Type_create_subarray( 1, headerSize, headerSize, headerStart, MPI_ORDER_C, MPI_INT, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Commit file type for reading file header
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// Open file
		fileOpenError = MPI_File_open( MPI_COMM_WORLD, (char*)fileName,
				                  	   MPI_MODE_RDONLY,
				                  	   MPI_INFO_NULL, &dataFile );

		// Set file view for reading header
		viewSetError = MPI_File_set_view( dataFile, 0, MPI_INT, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Read file header
		fileReadError = MPI_File_read_at_all( dataFile, 0, dataDim, nDim, MPI_INT, &status );
		assert( fileReadError == MPI_SUCCESS );

		// Free file type for header
		MPI_Type_free( &fileType );

		// Convert 1D myId to 3D block ID
		int quotient = myId;
		for( int i=0; i<nDim; i++ )
		{
			blockId[i] = quotient % nBlocks[i];
			quotient =  quotient / nBlocks[i];
		}

		// Compute size and bounding box of each block
		assert( ITL_util<int>::prod( nBlocks, nDim ) == nProcs );
		for( int i=0; i<nDim; i++ )
		{
			//blockDim[i] = (int) floor( (float)dataDim[i] / nBlocks[i] );
			blockDim[i] = (int) ceil( (float)dataDim[i] / nBlocks[i] );
			low[i] = blockId[i] * blockDim[i];

			// Make necessary adjustments if the last processor gets a block
			// of different size
			if( dataDim[i] % nBlocks[i] != 0 && blockId[i] == nBlocks[i]-1 )
				//blockDim[i] = blockDim[i] + dataDim[i] % nBlocks[i];
				blockDim[i] = dataDim[i] -low[i];

			high[i] = ( low[i] + blockDim[i]-1 );
			paddedLow[i] = ITL_util<int>::clamp( low[i] - neighborhoodSize, 0, dataDim[i]-1 );
			paddedHigh[i] = ITL_util<int>::clamp( high[i] + neighborhoodSize, 0, dataDim[i]-1 );
			paddedBlockDim[i] = paddedHigh[i] - paddedLow[i] + 1;

			lowPad[i] = low[i] - paddedLow[i];
			highPad[i] = paddedHigh[i] - high[i];

		}// end for

		// Compute number of elements in data and in this block
		int nel = ITL_util<int>::prod( dataDim, nDim );
		int nelBlock = ITL_util<int>::prod( blockDim, nDim );
		int nelBlockWithPad = ITL_util<int>::prod( paddedBlockDim, nDim );

		#ifdef DEBUG_MODE
		if( nDim == 3 )
		{
			printf( "%d: NX: %d, NY: %d, NZ: %d\nTotal number of elements in data: %d\n", myId, dataDim[0], dataDim[1], dataDim[2], nel );
			printf( "%d: Block ID: %d %d %d\n", myId, blockId[0], blockId[1], blockId[2] );
			printf( "%d: Block Start point: %d %d %d\n", myId, low[0], low[1], low[2] );
			printf( "%d: Block End point: %d %d %d\n", myId, high[0], high[1], high[2] );
			printf( "%d: Block size: %d, %d, %d\n", myId, blockDim[0], blockDim[1], blockDim[2] );
			printf( "%d: Block Start point (with pad): %d %d %d\n", myId, paddedLow[0], paddedLow[1], paddedLow[2] );
			printf( "%d: Block End point (with pad): %d %d %d\n", myId, paddedHigh[0], paddedHigh[1], paddedHigh[2] );
			printf( "%d: Block Size (with pad): %d %d %d\n", myId, paddedBlockDim[0], paddedBlockDim[1], paddedBlockDim[2] );
			printf( "%d: Lower side Pad Size: %d %d %d\n", myId, lowPad[0], lowPad[1], lowPad[2] );
			printf( "%d: Higher side Pad Size: %d %d %d\n", myId, highPad[0], highPad[1], highPad[2] );
			printf( "%d: Block id: %d, %d, %d\n", myId, blockId[0], blockId[1], blockId[2] );
		}
		#endif

		// Allocate memory for data to be read
		T* array = new T[nelBlockWithPad];

		// Header size
		MPI_Offset headerOffset = nDim * sizeof(MPI_INT);

		// Create eType for reading data
		MPI_Type_contiguous( sizeof(T), MPI_BYTE, &eType );
		commitError = MPI_Type_commit( &eType );
		assert( commitError == MPI_SUCCESS );

		// Specify sub-array to be written by this processor
		subarrayError = MPI_Type_create_subarray( nDim, dataDim, paddedBlockDim, paddedLow, MPI_ORDER_FORTRAN, eType, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Create file type for reading data
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// A few basic checks
		assert( array != NULL );
		assert( dataFile != NULL );
		assert( fileOpenError == MPI_SUCCESS );

		// Create file view
		viewSetError = MPI_File_set_view( dataFile, headerOffset, eType, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Write data to the current file view
		fileReadError = MPI_File_read_at_all( dataFile, 0, array, nelBlockWithPad, eType, &status );
		assert( fileReadError == MPI_SUCCESS );

		// close file
		MPI_File_close( &dataFile );

		// Free
		delete [] headerSize;
		delete [] headerStart;

		return array;

	}// end function


	/**
	 * Modified version of parallel file writer.
	 * Stores field as binary file in vec format. Uses advanced MPI-2 features like
	 * subarray, file view etc.
	 * @param data Pointer to portion of data
	 * @param fileName File name.
	 * @param dataDim Length of field along each dimension (can be a zero filled array).
	 * @param blockDim Length of the block along each dimension to be read by this instance
	 * of the program (can be a zero filled array).
	 * @param nBlocks Pointer to array containing number of blocks along each dimension.
	 * @param blockId Pointer to array containing serial number of block along each dimension.
	 * @param low Pointer to array containing lower grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param nDim Dimensionality of the field.
	 * @param myId Rank of this processor
	 * @param nProcs Total number of processes
	 */
	static void writeFieldBinaryParallel2( T* data, const char* fileName, int* dataDim,
									int* blockDim, int* nBlocks, int* blockId,
									float* low, int nDim,
									int myId, int nProcs )
	{
		MPI_File dataFile;
		MPI_Datatype eType, fileType;
		int fileOpenError, fileWriteError, subarrayError, viewSetError, commitError;
		MPI_Status status;
		MPI_Offset offset;

		// Compute total number of data points and the number of data points
		// to be written by this processor
		int nel = ITL_util<int>::prod( dataDim, nDim );
		int nelBlock = ITL_util<int>::prod( blockDim, nDim );

		int* lowInt = new int[nDim];
		for( int i=0; i<nDim; i++ )
			lowInt[i] = (int)floor( low[i] );

		#ifdef DEBUG_MODE
		if( blockDim[0] == 0 || blockDim[1] == 0 || blockDim[2] == 0 )
		{
			printf( "%d: Block size: %d, %d, %d\n", myId, blockDim[0], blockDim[1], blockDim[2] );
			printf( "%d: Block low limit: %d, %d, %d\n", myId, lowInt[0], lowInt[1], lowInt[2] );
			printf( "%d: Data size: %d, %d, %d\n", myId, dataDim[0], dataDim[1], dataDim[2] );
			printf( "%d: Block id: %d, %d, %d\n", myId, blockId[0], blockId[1], blockId[2] );
		}
		#endif

		// Specify sub-array for writing header
		int *headerSize = new int[1]; headerSize[0] = 3;
		int *headerStart = new int[1]; headerStart[0] = 0;
		subarrayError = MPI_Type_create_subarray( 1, headerSize, headerSize, headerStart, MPI_ORDER_C, MPI_INT, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Commit new file type
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// Open file
		fileOpenError = MPI_File_open( MPI_COMM_WORLD, (char*)fileName,
				                  	   MPI_MODE_WRONLY | MPI_MODE_CREATE,
				                  	   MPI_INFO_NULL, &dataFile );


		// Set file view for writing header
		viewSetError = MPI_File_set_view( dataFile, 0, MPI_INT, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Write data to the current file view
		fileWriteError = MPI_File_write_at_all( dataFile, 0, dataDim, nDim, MPI_INT, &status );
		assert( fileWriteError == MPI_SUCCESS );

		// Free file type
		MPI_Type_free( &fileType );

		// Header size
		MPI_Offset headerOffset = nDim * sizeof(MPI_INT);

		// Create eType
		MPI_Type_contiguous( sizeof(T), MPI_BYTE, &eType );
		commitError = MPI_Type_commit( &eType );
		assert( commitError == MPI_SUCCESS );

		// Specify sub-array for writing data
		subarrayError = MPI_Type_create_subarray( nDim, dataDim, blockDim, lowInt, MPI_ORDER_FORTRAN, eType, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Create file type for writing data
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// A few basic checks
		assert( data != NULL );
		assert( dataFile != NULL );
		assert( fileOpenError == MPI_SUCCESS );

		// Create file view for writing data
		viewSetError = MPI_File_set_view( dataFile, headerOffset, eType, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Write data to the current file view
		fileWriteError = MPI_File_write_at_all( dataFile, 0, data, nelBlock, eType, &status );
		assert( fileWriteError == MPI_SUCCESS );

		// close file
		MPI_File_close( &dataFile );

		// Free
		delete [] headerSize;
		delete [] headerStart;
		delete [] lowInt;

	}// end function

	/**
	 * Modified version of parallel file reader (allows different neighborhood size along different dimension).
	 * Stores field as binary file in vec format. Uses advanced MPI-2 features like
	 * subarray, file view etc.
	 * @param fileName File name.
	 * @param nDim Number of dimensions of data.
	 * @param dataDim Length of field along each dimension (can be a zero filled array).
	 * @param blockDim Length of the block along each dimension to be read by this instance
	 * of the program (can be a zero filled array).
	 * @param nBlocks Pointer to array containing number of blocks along each dimension.
	 * @param blockId Pointer to array containing serial number of block along each dimension.
	 * @param low Pointer to array containing lower grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param high Pointer to array containing higher grid (associated to the field) bounds in
	 * continuous space along each dimension (can be a zero filled array).
	 * @param lowPad Pointer to array containing ghost layer span along each dimension on the
	 * lower end (can be a zero filled array).
	 * @param highPad Pointer to array containing ghost layer span along each dimension on the
	 * upper end (can be a zero filled array).
	 * @param neighborhoodSize Pointer to neighborhood length along different directions for each point.
	 * @param myId Rank of this processor
	 * @param nProcs Total number of processes
	 * @return pointer to portion of data
	 */
	static T* readFieldBinaryParallel3( const char* fileName, int nDim, int* dataDim,
								int* blockDim,  int* nBlocks, int* blockId,
								int* low, int* high,
								int* lowPad, int* highPad, int* neighborhoodSize,
								int myId, int nProcs )
	{
		MPI_File dataFile;
		MPI_Datatype eType, fileType;
		int fileOpenError, fileReadError, subarrayError, viewSetError, commitError;
		MPI_Status status;
		MPI_Offset offset;
		int* paddedLow = new int[nDim];
		int* paddedHigh = new int[nDim];
		int* paddedBlockDim = new int[nDim];

		// Create subarray for file header
		int *headerSize = new int[1];  headerSize[0] = 3;
		int *headerStart = new int[1]; headerStart[0] = 0;
		subarrayError = MPI_Type_create_subarray( 1, headerSize, headerSize, headerStart, MPI_ORDER_C, MPI_INT, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Commit file type for reading file header
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// Open file
		fileOpenError = MPI_File_open( MPI_COMM_WORLD, (char*)fileName,
				                  	   MPI_MODE_RDONLY,
				                  	   MPI_INFO_NULL, &dataFile );

		// Set file view for reading header
		viewSetError = MPI_File_set_view( dataFile, 0, MPI_INT, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Read file header
		fileReadError = MPI_File_read_at_all( dataFile, 0, dataDim, nDim, MPI_INT, &status );
		assert( fileReadError == MPI_SUCCESS );

		// Free file type for header
		MPI_Type_free( &fileType );

		// Convert 1D myId to 3D block ID
		int quotient = myId;
		for( int i=0; i<nDim; i++ )
		{
			blockId[i] = quotient % nBlocks[i];
			quotient =  quotient / nBlocks[i];
		}

		// Compute size and bounding box of each block
		assert( ITL_util<int>::prod( nBlocks, nDim ) == nProcs );
		for( int i=0; i<nDim; i++ )
		{
			//blockDim[i] = (int) floor( (float)dataDim[i] / nBlocks[i] );
			blockDim[i] = (int) ceil( (float)dataDim[i] / nBlocks[i] );
			low[i] = blockId[i] * blockDim[i];

			// Make necessary adjustments if the last processor gets a block
			// of different size
			if( dataDim[i] % nBlocks[i] != 0 && blockId[i] == nBlocks[i]-1 )
				//blockDim[i] = blockDim[i] + dataDim[i] % nBlocks[i];
				blockDim[i] = dataDim[i] -low[i];

			high[i] = ( low[i] + blockDim[i]-1 );
			paddedLow[i] = ITL_util<int>::clamp( low[i] - neighborhoodSize[i], 0, dataDim[i]-1 );
			paddedHigh[i] = ITL_util<int>::clamp( high[i] + neighborhoodSize[i], 0, dataDim[i]-1 );
			paddedBlockDim[i] = paddedHigh[i] - paddedLow[i] + 1;

			lowPad[i] = low[i] - paddedLow[i];
			highPad[i] = paddedHigh[i] - high[i];

		}// end for

		// Compute number of elements in data and in this block
		int nel = ITL_util<int>::prod( dataDim, nDim );
		int nelBlock = ITL_util<int>::prod( blockDim, nDim );
		int nelBlockWithPad = ITL_util<int>::prod( paddedBlockDim, nDim );

		#ifdef DEBUG_MODE
		if( nDim == 3 )
		{
			//printf( "%d: NX: %d, NY: %d, NZ: %d\nTotal number of elements in data: %d\n", myId, dataDim[0], dataDim[1], dataDim[2], nel );
			//printf( "%d: Block ID: %d %d %d\n", myId, blockId[0], blockId[1], blockId[2] );
			//printf( "%d: Block Start point: %d %d %d\n", myId, low[0], low[1], low[2] );
			//printf( "%d: Block End point: %d %d %d\n", myId, high[0], high[1], high[2] );
			//printf( "%d: Block size: %d, %d, %d\n", myId, blockDim[0], blockDim[1], blockDim[2] );
			//printf( "%d: Block Start point (with pad): %d %d %d\n", myId, paddedLow[0], paddedLow[1], paddedLow[2] );
			//printf( "%d: Block End point (with pad): %d %d %d\n", myId, paddedHigh[0], paddedHigh[1], paddedHigh[2] );
			printf( "%d: Block Size (with pad): %d %d %d\n", myId, paddedBlockDim[0], paddedBlockDim[1], paddedBlockDim[2] );
			//printf( "%d: Lower side Pad Size: %d %d %d\n", myId, lowPad[0], lowPad[1], lowPad[2] );
			//printf( "%d: Higher side Pad Size: %d %d %d\n", myId, highPad[0], highPad[1], highPad[2] );
			//printf( "%d: Block id: %d, %d, %d\n", myId, blockId[0], blockId[1], blockId[2] );
		}
		#endif

		// Allocate memory for data to be read
		T* array = new T[nelBlockWithPad];

		// Header size
		MPI_Offset headerOffset = nDim * sizeof(MPI_INT);

		// Create eType for reading data
		MPI_Type_contiguous( sizeof(T), MPI_BYTE, &eType );
		commitError = MPI_Type_commit( &eType );
		assert( commitError == MPI_SUCCESS );

		// Specify sub-array to be written by this processor
		subarrayError = MPI_Type_create_subarray( nDim, dataDim, paddedBlockDim, paddedLow, MPI_ORDER_FORTRAN, eType, &fileType );
		assert( subarrayError == MPI_SUCCESS );

		// Create file type for reading data
		commitError = MPI_Type_commit( &fileType );
		assert( commitError == MPI_SUCCESS );

		// A few basic checks
		assert( array != NULL );
		assert( dataFile != NULL );
		assert( fileOpenError == MPI_SUCCESS );

		// Create file view
		viewSetError = MPI_File_set_view( dataFile, headerOffset, eType, fileType, "native", MPI_INFO_NULL );
		assert( viewSetError == MPI_SUCCESS );

		// Write data to the current file view
		fileReadError = MPI_File_read_at_all( dataFile, 0, array, nelBlockWithPad, eType, &status );
		assert( fileReadError == MPI_SUCCESS );

		// close file
		MPI_File_close( &dataFile );

		// Free
		delete [] headerSize;
		delete [] headerStart;

		return array;

	}// end function

	static void truncateFileHeader( char* inFileName, char* outFileName, int nDim )
	{
		// Open both files
		FILE* inFile = fopen( inFileName, "rb" );
		FILE* outFile = fopen( outFileName, "wb" );
		assert( inFile != NULL );
		assert( outFile != NULL );

		// Read header from input file
		char buffer[4];
		int dim[nDim];
		for( int i=0; i<nDim; i++ )
		{
			fread( buffer, 1, 4, inFile );
			dim[i] = *((int*)buffer);
		}
		// Compute total number of elements
		int nel = ITL_util<int>::prod( dim, nDim );

		// Read entire data to an array and write to the output file at once
		char *dataBuffer = new char[nel*sizeof(T)];
		fread( dataBuffer, sizeof(T), nel, inFile );
		fwrite( dataBuffer, sizeof(T), nel, outFile );

		#ifdef DEBUG_MODE
		if( nDim == 3 )
			printf( "NX: %d, NY: %d, NZ: %d\nTotal number of elements: %d\n", dim[0], dim[1], dim[2], nel );
		printf( "%d values copied to the output file\n", nel );
		#endif

		// close file
		fclose( inFile );
		fclose( outFile );

		// free
		delete [] dataBuffer;
		dataBuffer = NULL;

	}// end function

};

#endif
/* ITL_IOUTIL_H_ */
