/*
 * ITL_histogrammapper.h
 *
 *  Created on: Apr 8, 2012
 *  Author: abon
 */

#ifndef ITL_HISTOGRAMMAPPER_H_
#define ITL_HISTOGRAMMAPPER_H_

#include "ITL_util.h"
#include "ITL_histogram.h"
#include "ITL_trianglepatch.h"
//#include "ITL_geodesictree.h"
#include "ITL_field_regular.h"

template <class T>
class ITL_histogrammapper
{
	ITL_histogram* histogram;
	bool histogramRangeSet;
	SCALAR histogramMin, histogramMax;
	SCALAR histogramMinArray[2], histogramMaxArray[2];
	list<VECTOR3> vertexList[20];
	VECTOR3** vertexList2;
	list<ITL_trianglepatch> triangleList[20];

	float* theta;					/**< Angle variable. Angle related to spherical coordinates. */
	float* phi;						/**< Angle variable. Angle related to spherical coordinates. */

	bool bIsAngleMapInitialized;    /**< Boolean flag. Flag is set if the angle maps have been initialized. */
	bool bIsPatchRead;				/**< Boolean flag. Flag is set if the patch file has been read. */

	int iNrOfThetas;				/**< Histogram parameter. Number of thetas.  */
	int iNrOfPhis;					/**< Histogram parameter. Number of phis.  */
	float fNrOfThetas;				/**< Histogram parameter. Number of phis.  */
	float fNrOfPhis;				/**< Histogram parameter. Number of phis.  */
	int* piAngleMap;				/**< Histogram parameter. Mapping from vector to a region in sperical coordinates.  */



public:

	ITL_histogrammapper( ITL_histogram* hist )
	{
		this->histogram = hist;
		histogramRangeSet = false;
	}

	void
	setHistogramRange( T m, T M )
	{
		histogramMin = m;
		histogramMax = M;
		histogramRangeSet = true;
	}

	/**
	 * Joint histogram range set function.
	 * Sets the range for a 2D joint histogram.
	 * @param m1 Lower limit of the range of the fist histogram
	 * @param M1 Upper limit of the range of the fist histogram
	 * @param m2 Lower limit of the range of the second histogram
	 * @param m2 Upper limit of the range of the second histogram
	 */
	void
	setJointHistogramRange( T m1, T M1, T m2, T M2 )
	{
		histogramMinArray[0] = m1;
		histogramMaxArray[0] = M1;
		histogramMinArray[1] = m2;
		histogramMaxArray[1] = M2;
		histogramRangeSet = true;
	}

	/**
	 * Histogram bin assignment function for scalar fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Scalar( ITL_field_regular<T>* dataField,
									 ITL_field_regular<int>** binField,
									 int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		SCALAR nextV, rangeValue;

		//cout << "in" << endl;
		assert( dataField->getDataFull() != NULL );
		//cout << "in2" << endl;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}
		//cout << "in3" << endl;

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			//cout << "setting histogram range" << endl;
			// Get min-max values of the scalar field
			// ******************* Need to resolve typecast issue begin *************************
			//minValue = ITL_util<SCALAR>::Min( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			//maxValue = ITL_util<SCALAR>::Max( (SCALAR*)this->dataField->datastore->array, this->dataField->grid->nVertices );
			histogramMin = ITL_util<SCALAR>::Min( (SCALAR*)dataField->getDataFull(), dataField->getSize() );
			histogramMax = ITL_util<SCALAR>::Max( (SCALAR*)dataField->getDataFull(), dataField->getSize() );

			// ******************* Need to resolve typecast issue end *************************
		}
		rangeValue = histogramMax - histogramMin;

		// Compute bin width
		float binWidth = rangeValue / (float)nBin;

		#ifdef DEBUG_MODE
		printf( "Min: %g Max: %g Range: %g of the scalar values\n", histogramMin, histogramMax, rangeValue );
		printf( "Binwidth: %g\n", binWidth );
		#endif

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		int dimWithPad[4];
		(*binField)->getSizeWithPad( dimWithPad );
		//cout << "starting iteration" << endl;
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get scalar value at location
					// ******************* Need to resolve typecast issue *************************
					nextV = dataField->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV - histogramMin ) / binWidth  );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binField)->setDataAt( index1d, binId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}
		///cout << "data scan completed" << endl;

	}// end function


	/**
	 * Histogram bin assignment function for vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Vector( ITL_field_regular<T>* dataField,
			 	 	 	 	 	 	 ITL_field_regular<int>** binField,
			 	 	 	 	 	 	 int nBin, int iRes = 0, char* binMapFile = NULL )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField->getDataFull() != NULL );
		VECTOR3 nextV;

		//cout << "in" << endl;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}
		//cout << "in3" << endl;

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		(*binField)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", dimWithPad[0], dimWithPad[1], dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					nextV = (VECTOR3)dataField->getDataAt( index1d );
					//cout << nextV[0] << " " << nextV[1] << " " << nextV[2] << endl;

					// Obtain the binID corresponding to the value at this location
					binId = histogram->get_bin_number_3D( nextV, iRes );
					//cout << binId << endl;
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					//cout << binId << endl;
					(*binField)->setDataAt( index1d, binId );



					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}
		//cout << "data scan completed" << endl;
        // delete lPadHisto;
        // delete hPadHisto;

	}// end function

	// Read the lookup table provided in the header file
	void
	readPatches_header( int nBin )
	{
		static char szPatch[] = {
		#include "ITL_patch.h"
		};

		char *szToken = strtok(szPatch, "\n");
		float f2Temp[2];

		for(int i=0;i<nBin; i++)
		{
			float fAngle_radian;
			sscanf( szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
			theta[i*2+0] = f2Temp[0];
			theta[i*2+1] = f2Temp[1];
			szToken = strtok(NULL, "\n");
			//printf( "%f %f\n", f2Temp[0], f2Temp[1] );

			sscanf(szToken,"%f, %f", &f2Temp[0],&f2Temp[1] );
			phi[i*2+0] = f2Temp[0];
			phi[i*2+1] = f2Temp[1];
			szToken = strtok(NULL, "\n");
			//printf( "%f %f\n", f2Temp[0], f2Temp[1] );
		}
	}

	/**
	 * Histogram bin assignment function for vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Vector_Geodesic( ITL_field_regular<T>* dataField,
			 	 	 	 	 	 	 	      ITL_field_regular<int>** binField,
			 	 	 	 	 	 	 	      int nBin, int nDivision )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField->getDataFull() != NULL );
		VECTOR3 nextV;

		//cout << "in" << endl;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}
		//cout << "in3" << endl;

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		(*binField)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", dimWithPad[0], dimWithPad[1], dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					nextV = (VECTOR3)dataField->getDataAt( index1d );
					//cout << nextV[0] << " " << nextV[1] << " " << nextV[2] << endl;

					// Obtain the binID corresponding to the value at this location
					//binId = getBinNumber3D( nextV, &vertexList[nDivision-1], &triangleList[nDivision-1]  );
					binId = getBinNumber3DViaTable( nextV, nBin );
					//cout << binId << endl;
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					//cout << binId << endl;
					(*binField)->setDataAt( index1d, binId );



					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}
		//cout << "data scan completed" << endl;
        // delete lPadHisto;
        // delete hPadHisto;

	}// end function

	/**
	 * Histogram bin assignment function for 2D vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Vector2( ITL_field_regular<T>* dataField,
	   	 	 	 	 	 	 	      ITL_field_regular<int>** binField,
	 	 	 	 	 	 	 	 	  int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField->getDataFull() != NULL );
		T *nextV = new T();

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		(*binField)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", this->dataField->grid->dimWithPad[0], this->dataField->grid->dimWithPad[1], this->dataField->grid->dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					*nextV = dataField->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location
					binId = histogram->get_bin_number_2D( *nextV, nBin );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binField)->setDataAt( index1d, binId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

        delete nextV;

	}// end function

	/**
	 * Histogram bin assignment function for vector fields.
	 * Creates a scalar field of histogram at each grid vertex.
	 * @param nBins Number of bins to use in histogram computation.
	 */
	void
	computeHistogramBinField_Vector_Decomposed( ITL_field_regular<T>* dataField,
			 	 	 	 	 	 	 	 	 	ITL_field_regular<int>** binFieldX,
			 	 	 	 	 	 	 	 	 	ITL_field_regular<int>** binFieldY,
			 	 	 	 	 	 	 	 	 	ITL_field_regular<int>** binFieldZ,
			 	 	 	 	 	 	 	 	 	float* mMArray,
			 	 	 	 	 	 	 	 	 	int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		// Compute bin width
		float binWidthX = ( mMArray[1] - mMArray[0] ) / (float)nBin;
		float binWidthY = ( mMArray[3] - mMArray[2] ) / (float)nBin;
		float binWidthZ = ( mMArray[5] - mMArray[4] ) / (float)nBin;

		assert( dataField->getDataFull() != NULL );
		VECTOR3 nextV;

		// Initialize the padded scalar field for histogram bins
		if( (*binFieldX) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield-X Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield-X High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield-X Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield-X High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binFieldX) = new ITL_field_regular<int>( dataField->getNumDim(),
													   low, high,
													   lowPad, highPad,
													   neighborhoodSize );
		}
		if( (*binFieldY) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield-Y Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield-Y High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield-Y Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield-Y High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binFieldY) = new ITL_field_regular<int>( dataField->getNumDim(),
													   low, high,
													   lowPad, highPad,
													   neighborhoodSize );
		}
		if( (*binFieldZ) == NULL )
		{
			//cout << "creating new binfield" << endl;
			dataField->getBounds( low, high );
			dataField->getPadSize( lowPad, highPad );
			dataField->getNeighborhoodSize( neighborhoodSize );

			#ifdef DEBUG_MODE
			printf( "binfield-Z Low: %g %g %g\n", low[0], low[1], low[2] );
			printf( "binfield-Z High: %g %g %g\n", high[0], high[1], high[2] );
			printf( "binfield-Z Low Pad: %d %d %d\n", lowPad[0], lowPad[1], lowPad[2] );
			printf( "binfield-Z High Pad: %d %d %d\n", highPad[0], highPad[1], highPad[2] );
			#endif

			(*binFieldZ) = new ITL_field_regular<int>( dataField->getNumDim(),
													   low, high,
													   lowPad, highPad,
													   neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		(*binFieldX)->getSizeWithPad( dimWithPad );

		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					//printf( "%d %d %d %d %d %d %d\n", dimWithPad[0], dimWithPad[1], dimWithPad[2], x, y, z, index1d );
					// Get vector at location
					nextV = (VECTOR3)dataField->getDataAt( index1d );
					//nextV.Normalize();
					//cout << nextV[0] << " " << nextV[1] << " " << nextV[2] << endl;

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV.x() - mMArray[0] ) / binWidthX );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binFieldX)->setDataAt( index1d, binId );

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV.y() - mMArray[2] ) / binWidthY );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binFieldY)->setDataAt( index1d, binId );

					// Obtain the binID corresponding to the value at this location
					binId = (int)floor( ( nextV.z() - mMArray[4] ) / binWidthZ );
					binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
					(*binFieldZ)->setDataAt( index1d, binId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}
		//cout << "data scan completed" << endl;
        // delete lPadHisto;
        // delete hPadHisto;

	}// end function


	void
	computeJointHistogramBinField_Scalar( ITL_field_regular<T>* dataField1,
	  								      ITL_field_regular<T>* dataField2,
	 	 	 	   	   	   	   	   	   	  ITL_field_regular<int>** binField,
	 	 	 	   	   	   	   	   	   	  int nBin )
	{
		cout << "Scalar" << endl;
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		T nextVField1, nextVField2;

		//assert( this->dataField1->datastore->array != NULL );
		//assert( this->dataField2->datastore->array != NULL );
		assert( dataField1->getDataFull() != NULL );
		assert( dataField2->getDataFull() != NULL );

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField1->getBounds( low, high );
			dataField1->getPadSize( lowPad, highPad );
			dataField1->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField1->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Compute the range over which histogram computation needs to be done
		if( histogramRangeSet == false )
		{
			// Get min-max values of both the scalar fields
			histogramMinArray[0] = ITL_util<T>::Min( dataField1->getDataFull(), dataField1->getSize() );
			histogramMaxArray[0] = ITL_util<T>::Max( dataField1->getDataFull(), dataField1->getSize() );
			histogramMinArray[1] = ITL_util<T>::Min( dataField2->getDataFull(), dataField2->getSize() );
			histogramMaxArray[1] = ITL_util<T>::Max( dataField2->getDataFull(), dataField2->getSize() );
		}

		// Compute bin widths along each dimension
		T rangeValue1 = histogramMaxArray[0] - histogramMinArray[0];
		float binWidth1 = rangeValue1 / (float)nBin;
		T rangeValue2 = histogramMaxArray[1] - histogramMinArray[1];
		float binWidth2 = rangeValue2 / (float)nBin;

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId1, binId2, combinedBinId;
		int dimWithPad[4];
		dataField1->getSizeWithPad( dimWithPad );
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					nextVField1 = dataField1->getDataAt( index1d );
					nextVField2 = dataField1->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location from both fields individually
					binId1 = (int)floor( ( nextVField1 - histogramMinArray[0] ) / binWidth1  );
					binId1 = ITL_util<int>::clamp( binId1, 0, nBin-1 );
					binId2 = (int)floor( ( nextVField2 - histogramMinArray[1] ) / binWidth2  );
					binId2 = ITL_util<int>::clamp( binId2, 0, nBin-1 );

					// Combine the two bin indices from the two fields
					combinedBinId = binId2 * nBin + binId1;

					// Store the bin ID
					(*binField)->setDataAt( index1d, combinedBinId );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

		// delete lPadHisto;
		// delete hPadHisto;

	}// end function

	/**
	 * Joint histogram bin assignment function for vector field.
	 * @param nBin Number of bins to be used in joint histogram computation.
	 */
	void
	computeJointHistogramBinField_Vector( ITL_field_regular<T>* dataField1,
										  ITL_field_regular<T>* dataField2,
										  ITL_field_regular<int>** binField,
										  int nBin )
	{
		float low[4];
		float high[4];
		int lowPad[4];
		int highPad[4];
		int neighborhoodSize[4];
		int dimWithPad[4];

		assert( dataField1->getDataFull() != NULL );
		assert( dataField2->getDataFull() != NULL );
		VECTOR3 nextVField1;
		VECTOR3 nextVField2;

		// Initialize the padded scalar field for histogram bins
		if( (*binField) == NULL )
		{
			dataField1->getBounds( low, high );
			dataField1->getPadSize( lowPad, highPad );
			dataField1->getNeighborhoodSize( neighborhoodSize );

			(*binField) = new ITL_field_regular<int>( dataField1->getNumDim(),
													  low, high,
													  lowPad, highPad,
													  neighborhoodSize );
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		dataField1->getSizeWithPad( dimWithPad );
		for( int z=0; z<dimWithPad[2]; z++ )
		{
			for( int y=0; y<dimWithPad[1]; y++ )
			{
				for( int x=0; x<dimWithPad[0]; x++ )
				{
					// Get vector at this location from both fields individually
					nextVField1 = dataField1->getDataAt( index1d );
					nextVField2 = dataField2->getDataAt( index1d );

					// Obtain the binID corresponding to the value at this location from both fields individually
					int binValue1 = histogram->get_bin_number_3D( nextVField1 );
					int binValue2 = histogram->get_bin_number_3D( nextVField2 );

					// Combine the two bin indices from the two fields
					int combinedBin = binValue2 * nBin + binValue1;

					// Store the bin ID
					(*binField)->setDataAt( index1d, combinedBin );

					// increment to the next grid vertex
					index1d += 1;
				}
			}
		}

	}// end function

	static void
	computeHistogramFrequencies( ITL_field_regular<int>** binField,
								 int* freqList,
								 int nBin )
	{
		assert( binField != NULL );
		assert( freqList != NULL );

		// Scan through bin Ids and keep count
		for( int i=0; i<nBin; i++ )
			freqList[i] = 0;

		for( int i=0; i<(*binField)->getSize(); i++ )
		{
			//cout << (*binField)->getDataAt(i) << endl;
			freqList[ (*binField)->getDataAt(i) ] ++;
		}
	}// end function

	static void
	computeHistogramFrequencies( ITL_field_regular<int>** binField,
								 float* freqList,
								 int nBin )
	{
		assert( binField != NULL );
		assert( freqList != NULL );

		// Scan through bin Ids and keep count
		for( int i=0; i<nBin; i++ )
			freqList[i] = 0;

		for( int i=0; i<(*binField)->getSize(); i++ )
		{
			//cout << (*binField)->getDataAt(i) << endl;
			freqList[ (*binField)->getDataAt(i) ] ++;
		}

		for( int i=0; i<nBin; i++ )
			freqList[i] /= (float)(*binField)->getSize();

	}// end function

	static void
	computeHistogramFrequencies( int *binIds, int* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
	}// end function

	static void
	computeHistogramFrequencies( int *binIds, float* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
	}// end function

	static void
	computeHistogramFrequencies( int *binIds, double* freqList, int nPoint )
	{
		// Scan through bin Ids and keep count
		for( int i=0; i<nPoint; i++ )
			freqList[ binIds[i] ] ++;
	}// end function

	/**
	  * Cross-validation function.
	  */
	void crossValidate_scalar( ITL_field_regular<SCALAR>* dataField,
			  				   int nIter, int start, int step,
			  				   int *nBinArray,
			  				   double *scoreList,
			  				   int* nBinOptimal )
	{
		int nDim = dataField->getNumDim();
		float lowF[3], highF[3];
		dataField->getBounds( lowF, highF );
		double N = dataField->getSize();
		ITL_field_regular<int>* binField = NULL;
		int *freqListPtr = NULL;
		double *normFreqListPtr = NULL;

		// Iterate
		double minScore = 10000000;
		int minScoreIndex = -1;
		double h = 0;

		// Scan through different number of bins
		for( int i = 0; i<nIter; i++ )
		{
			// Compute number of bins
			nBinArray[i] = start + i*step;

			// Create bin field
			binField = new ITL_field_regular<int>( nDim, lowF, highF );

			// Compute bin field
			computeHistogramBinField_Scalar( dataField, &binField, nBinArray[i] );

			// Compute frequencies
			freqListPtr = new int[nBinArray[i]];
			memset( freqListPtr, 0, sizeof(int)*nBinArray[i] );
			computeHistogramFrequencies( binField->getDataFull(),
					 	 	 	 	 	 freqListPtr,
										 binField->getSize() );

			// Compute bin width
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			// Normalize frequencies
			normFreqListPtr = new double[nBinArray[i]];
			memset( normFreqListPtr, 0, sizeof(double)*nBinArray[i] );
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr[j] / N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			//#endif

			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			delete [] normFreqListPtr;
			delete [] freqListPtr;
			delete binField;

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];

	}// End function


	/** added by Tzu-Hsuan
	 *
	 */
	void
	crossValidate_optimized_scalar( ITL_field_regular<T>* dataField,
	 					  	  	    int nIter, int nBinMax,
	 					  	  	    int *nBinArray,
	 					  	  	    double* scoreList,
	 					  	  	    int* nBinOptimal )
	{
		int nDim = dataField->getNumDim();
		float lowF[3], highF[3];
		dataField->getBounds( lowF, highF );
		double N = dataField->getSize();
		//printf( "N: %f\n", N );
		int freqListPtr[nBinMax];
		int freqListPtr2[nBinMax];
		double normFreqListPtr[nBinMax];
		int cross_nBin;
		double minScore;
		int minScoreIndex;
		double h = 0;

		// Open file for writing
		FILE* outFile = NULL;
		stringstream temp( stringstream::in | stringstream::out );

		// Create bin field
		ITL_field_regular<int>* binField = new ITL_field_regular<int>( nDim, lowF, highF );

		// Iterate over a range of number of bins
		for( int i = nIter-1; i>=0; i-- )
		{
			// Update number of bins
			if( i == nIter-1 )
				nBinArray[i] = nBinMax;
			else
				nBinArray[i] = nBinArray[i+1]/2;

			// compute/update frequencies
			if( i == nIter-1 )
			{
				// Compute bin field from data only in the first iteration
				computeHistogramBinField_Scalar( dataField, &binField, nBinArray[i] );

				// Compute frequencies
				memset( freqListPtr2, 0, sizeof(int)*nBinArray[i] );
				ITL_histogrammapper::computeHistogramFrequencies( binField->getDataFull(), freqListPtr2,  binField->getSize() );

			}
			else
			{
				memset( freqListPtr2, 0, sizeof(int)*nBinArray[i] );
				for( int j=0; j<=nBinArray[i]; j++ )
					freqListPtr2[j] = freqListPtr[j*2] + freqListPtr[j*2+1];
			}

			// Normalize updated frequencies
			memset( normFreqListPtr, 0, sizeof(double)*nBinArray[i] );
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr2[j] / N;

			// Update bin width
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			stringstream temp( stringstream::in | stringstream::out );
			temp << "/home/abon/histotest" << i << ".csv";
			string filename = temp.str();
			outFile = fopen( filename.c_str(), "w" );

			for( int k=0; k<nBinArray[i]; k++ )
				fprintf( outFile, "%d, %g\n", freqListPtr2[k], normFreqListPtr[k] );
			fclose( outFile );

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			//#endif

			// Update min score
			if( i == nIter-1 || scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			// Store the current frequent list for next iteration
			memcpy( freqListPtr, freqListPtr2, sizeof(int)*nBinArray[i] );

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];


	}// End function

	void
	crossValidate_optimized2( int* freqList,
							  int nIter, int nBinMax, int N,
							  int *nBinArray,
  							  double* scoreList,
  							  int* nBinOptimal )
	{
		int* freqListPtr = new int[nBinMax];
		memcpy( freqListPtr, freqList, sizeof(int)*nBinMax );
		int* freqListPtr2 = NULL;
		double* normFreqListPtr = NULL;
		double h = 0;

		// Compute number of bins for different iterations
		for( int i = 0; i<nIter; i++ )
			nBinArray[nIter-i-1] = (int)(nBinMax/(pow(2.0,(double)i)));

		// Normalize frequencies
		normFreqListPtr = new double[nBinArray[nIter-1]];
		for( int j=0; j<nBinArray[nIter-1]; j++ )
			normFreqListPtr[j] =  freqListPtr[j] / (double)N;

		// Compute cross-validation score
		h = ( histogramMax - histogramMin ) / (double) nBinArray[nIter-1];
		scoreList[nIter-1] = computeCrossValidationScore( normFreqListPtr, nBinArray[nIter-1], N, h );

		//#ifdef DEBUG_MODE
		printf( "%d, %g, %g\n", nBinArray[nIter-1], h, scoreList[nIter-1] );
		//#endif

		// Record latest minimum score
		double minScore = scoreList[nIter-1];
		int minScoreIndex = (nIter-1);

		// Clear temp resources
		delete [] normFreqListPtr;

		// Update
		for( int i = nIter-2; i>=0; i-- )
		{
			// Update bin width
			h = ( histogramMax - histogramMin ) / (double) nBinArray[i];

			// Update frequencies
			freqListPtr2 = new int[nBinArray[i]];
			for( int j=0; j<=nBinArray[i]; j++ )
				freqListPtr2[j] = freqListPtr[j*2] + freqListPtr[j*2+1];

			// Normalize updated frequencies
			normFreqListPtr = new double[nBinArray[i]];
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListPtr2[j] / (double)N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d, %g, %g\n", nBinArray[i], h, scoreList[i] );
			//#endif

			if( scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			// Store the current frequent list for next iteration
			memcpy( freqListPtr, freqListPtr2, sizeof(int)*nBinArray[i] );

			delete [] normFreqListPtr;
			delete [] freqListPtr2;

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];

		delete [] freqListPtr;

	}// End function

	void
	crossValidate_optimized_vector( int* freqList,
								    double N, int nBinMax,
								    int nDivision,
								    int* nBinArray, double* scoreList,
								    int* nBinOptimal )
	{
		//int nBinArray[nDivision];
		//double scoreList[nDivision];
		int* freqListCur = NULL;
		int* freqListLast = NULL;
		double* normFreqListPtr = NULL;
		double h = sqrt( (4*pi) / (double)nBinMax );
		double minScore = scoreList[nDivision-1];
		int minScoreIndex = (nDivision-1);
		list<ITL_trianglepatch>::iterator iterLast;
		list<ITL_trianglepatch>::iterator iterCur;

		// Iterate and update
		// We are moving from child to the parent levels
		// cur -- parent, last -- child level
		for( int i = nDivision-1; i>=0; i-- )
		{
			// Calculate number of bins
			if( i == nDivision-1 )
				nBinArray[i] = nBinMax;
			else
				nBinArray[i] = nBinArray[i+1] / 4;

			// Calculate bin width
			//h = pow( 4.0, i);
			h = h * 4.0;

			// Update frequencies
			if( i == nDivision-1 )
			{
				freqListCur = new int[nBinArray[i]];
				freqListLast = new int[nBinArray[i]];
				memcpy( freqListCur, freqList, sizeof( int )*nBinArray[i] );
				memcpy( freqListLast, freqList, sizeof( int )*nBinArray[i] );
			}
			else
			{
				delete [] freqListCur;
				freqListCur = new int[nBinArray[i]];
				memset( freqListCur, 0, sizeof( int )*nBinArray[i] );

				int iLast = 0;
				//int nv = triangleList[i+1].size();
				//printf( "Number of bins in child level: %d\n", nv );
				for( iterLast = triangleList[i+1].begin();
					 iterLast != triangleList[i+1].end();
					 iterLast++ )
				{
					int parentId = iterLast->getParentID();
					//cout << parentId << " " << freqListLast[iLast] << endl;
					freqListCur[parentId] += freqListLast[iLast];
					iLast++;
				}

				//for( int l=0; l<nBinArray[i]; l++ )
				//	printf( "%d, ", freqListCur[l] );
				//printf( "\n" );

				// Will be useful in the next iteration
				delete [] freqListLast;
				freqListLast = new int[nBinArray[i]];
				memcpy( freqListLast, freqListCur, sizeof( int )*nBinArray[i] );

			}

			// Normalize updated frequencies
			normFreqListPtr = new double[nBinArray[i]];
			memset( normFreqListPtr, 0, sizeof(double)*nBinArray[i] );
			for( int j=0; j<nBinArray[i]; j++ )
				normFreqListPtr[j] =  freqListCur[j] / N;

			// Compute cross-validation score
			scoreList[i] = computeCrossValidationScore( normFreqListPtr, nBinArray[i], N, h );

			//#ifdef DEBUG_MODE
			printf( "%d: %d, %g, %g\n", i, nBinArray[i], h, scoreList[i] );
			//#endif

			// Update minimum score
			if( i == nDivision-1 || scoreList[i] < minScore )
			{
				minScore = scoreList[i];
				minScoreIndex = i;
			}

			// Free temporary resources
			delete [] normFreqListPtr;

		}// end for : iterations

		(*nBinOptimal) = nBinArray[minScoreIndex];

	}// End function

	/**
	 * Cross-validation score computation function.
	 *
	 */
	double
	computeCrossValidationScore( double* fList, int nBin, double N, double h )
	{
		double sumOfSquares = 0;
		//double sum = 0;

		for( int j = 0; j<nBin; j++ )
		{
			assert( fList[j] >= 0 );
			sumOfSquares += ( fList[j] * fList[j] );
			//sum += fList[j];
			//if( nBin < 64 )
			//	printf("%.2f,", fList[j]);
		}
		//printf("\n");


		double cv_nom = 2 - sumOfSquares * (N+1.0);
		double cv_denom = (N-1.0) * h;
		double score = cv_nom / cv_denom;

		//printf( "d1 d2 d3: %g, %g, %g %g %g\n", N, h, cv_nom, cv_denom, score );

		//N = sum;
		//double firstPart =  2.0 / ((N - 1.0)*h);
		//double secondPart = (double)( N+1 ) / (double)(( N-1 )*h);

		//double denom1 = N*N;
		//double denom2 = N*N*h;
		//double denom3 = (double)(N+1)/(double)(N-1);
		//printf( "d1 d2 d3: %g %g %g\n", denom1, denom2, denom3 );
		//double secondPart = denom3 / denom2;

		//double secondPart = (double)( N+1 ) / (double)(( N-1 )*N*N*h);
		//double secondPart = (double)( N + 1.0 ) / (double)(( N - 1.0 )*h);
		//printf( "%g %g %g\n", firstPart, secondPart, sumOfSquares );
		//double score = firstPart - ( secondPart * sumOfSquares );

		return score;

	}// end function

	void
	createSphericalGrid( int nDivision, int* nBin )
	{
		assert( nDivision >= 1 );
		assert( nDivision <= 20 );
		(*nBin) = (int)( 20 * pow( 4.0, nDivision-1 ) );
		//printf( "Number of bins: %d\n", (*nBin) );

		// Calculate the golden ratio.
		double PHI = 2 * cos( pi / 5.0 );

		// The number of faces in an icosahedron.
		int NUM_TRIANGLES_IN_AN_ICOSAHEDRON = 20;

		// The dimensionality of the sphere.
		int NUM_DIMENSIONS = 3;

		// The number of vertices in a triangle.
		int NUM_VERTICES_IN_A_TRIANGLE = 3;

		// The number of triangles produced
		int NUM_NEW_TRIANGLES_PER_TRIANGLE = 4;

		// Vertices
		//VERTICES = [0, PHI, 1;
		//			  0, -PHI, 1;
		//			  0, PHI, -1;
		//			  0, -PHI, -1;
		//			  1, 0, PHI;
		//			  -1, 0, PHI;
		//			  1, 0, -PHI;
		//			  -1, 0, -PHI;
		//			  PHI, 1, 0;
		//			  -PHI, 1, 0;
		//			  PHI, -1, 0;
		//			  -PHI, -1, 0]
		//			  / norm([1, PHI]);

		double mg = sqrt( 1*1 + PHI*PHI );
		VECTOR3 v0( 0 / mg, PHI / mg, 1 / mg ); 	vertexList[0].push_back( v0 );
		VECTOR3 v1( 0 / mg, -PHI / mg, 1 / mg );	vertexList[0].push_back( v1 );
		VECTOR3 v2( 0 / mg, PHI / mg, -1 / mg );	vertexList[0].push_back( v2 );
		VECTOR3 v3( 0 / mg, -PHI / mg, -1 / mg );	vertexList[0].push_back( v3 );
		VECTOR3 v4( 1 / mg, 0 / mg, PHI / mg );		vertexList[0].push_back( v4 );
		VECTOR3 v5( -1 / mg, 0 / mg, PHI / mg );	vertexList[0].push_back( v5 );
		VECTOR3 v6( 1 / mg, 0 / mg, -PHI / mg );	vertexList[0].push_back( v6 );
		VECTOR3 v7( -1 / mg, 0 / mg, -PHI / mg );	vertexList[0].push_back( v7 );
		VECTOR3 v8( PHI / mg, 1 / mg, 0 / mg );		vertexList[0].push_back( v8 );
		VECTOR3 v9( -PHI / mg, 1 / mg, 0 / mg );	vertexList[0].push_back( v9 );
		VECTOR3 v10( PHI / mg, -1 / mg, 0 / mg );	vertexList[0].push_back( v10 );
		VECTOR3 v11( -PHI / mg, -1 / mg, 0 / mg );	vertexList[0].push_back( v11 );

		// Vertex indices
		int A[] = { 1, 4, 8, 6, 10, 3, 5, 1, 0, 4, 2, 8, 7, 6, 11, 11, 5, 0, 2, 9 };
		int B[] = { 3, 10, 4, 8, 6, 1, 11, 4, 5, 8, 0, 6, 2, 3, 7, 5, 0, 2, 9, 7 };
		int C[] = { 10, 1, 10, 10, 3, 11, 1, 5, 4, 0, 8, 2, 6, 7, 3, 9, 9, 9, 7, 11 };

		// Triangle indices
		for( int iT=0; iT<20; iT++ )
		{
			ITL_trianglepatch nextPatch( A[iT], B[iT], C[iT], 0 );
			triangleList[0].push_back( nextPatch );
		}

		// Iteratively create the triangle patches
		int level = 0;
		for( int i=0; i<nDivision-1; i++ )
		{
			// Split each triangle into 4 smaller triangles.
			divideTriangles2( &vertexList[i],
							 &triangleList[i],
							 &vertexList[i+1],
 							 &triangleList[i+1],
							 level,
							 NUM_DIMENSIONS,
							 NUM_VERTICES_IN_A_TRIANGLE,
							 NUM_NEW_TRIANGLES_PER_TRIANGLE );

			//cout << "iteration: " << i << endl;
			level ++;
		}

		vertexList2 = new VECTOR3*[nDivision];
		for( int i=0; i<nDivision; i++ )
		{
			// Copy vertices from list to array
			int nV = (int)vertexList[i].size();
			vertexList2[i] = new VECTOR3[nV];
			list<VECTOR3>::iterator iter = vertexList[i].begin();
			for( int iV=0; iV<nV; iV++ )
			{
				vertexList2[i][iV] = (*iter);
				iter++;
			}
		}

		// print (and save) Results
		FILE* outFile;
		int nv, nt;
		list<VECTOR3>::iterator iter1;
		list<ITL_trianglepatch>::iterator iter2;
		for( int i=0; i<nDivision; i++ )
		{
			// Open file for writing
			stringstream temp( stringstream::in | stringstream::out );
			temp << "/home/abon/spheregrid" << i << ".csv";
			string filename = temp.str();
			outFile = fopen( filename.c_str(), "w" );
			//cout << filename << endl;
			assert( outFile != NULL );

			nv = (int)vertexList[i].size();
			nt = (int)triangleList[i].size();
			//fprintf( outFile, "%d, %d, %d\n", nv, nt, 0 );

			//printf( "Level: %d\n", i );
			//printf( "Number of vertices: %d\n", nv );
			//printf( "Number of triangles: %d\n", nt );
			//printf( "" );

			iter1 = vertexList[i].begin();
			for( int j=0; j<nv; j++ )
			{
				fprintf( outFile, "%g, %g, %g\n", iter1->x(), iter1->y(), iter1->z() );
				iter1 ++;
			}

			iter2 = triangleList[i].begin();
			for( int j=0; j<nt; j++ )
			{
				fprintf( outFile, "%d, %d, %d\n", iter2->getVertexID(0), iter2->getVertexID(1), iter2->getVertexID(2) );
				iter2 ++;
			}

			// Close file
			fclose( outFile );
		}

		//FILE* files[nDivision];
		//string filename;
		//for( int i=0; i<nDivision; i++ )
		//{
			// Open file for writing
		//	stringstream temp( stringstream::in | stringstream::out );
		//	temp << "/home/abon/spheregrid" << i << ".csv";
		//	filename = temp.str();
		//	files[i] = fopen( filename.c_str(), "r" );
		//}

		FILE* writeFile = fopen( "/home/abon/spheregridhierarchy.csv", "w" );
		list<ITL_trianglepatch>::iterator iter3;
		int curId, pId;
		for( int i=0; i<nt; i++ )
		{
			// Save current triangle id
			fprintf( writeFile, "%d, ", i );
			curId = i;

			for( int j=nDivision-1; j>=1; j-- )
			{
				// Find parent ID of the current triangle
				iter3 = triangleList[j].begin();
				for( int k=0; k<curId; k++ )
					iter3++;
				pId = iter3->getParentID();

				// Store parent ID
				if( j > 1 )
					fprintf( writeFile, "%d, ", pId );
				else
					fprintf( writeFile, "%d\n", pId );

				// Parent ID becomes current for the next iteration
				curId = pId;
			}

		}

		fclose( writeFile );

		//for( int i=0; i<nDivision; i++ )
		//{
		//	fclose( files[i] );
		//}

	}// end function

	void
	divideTriangles2( list<VECTOR3>* lastVertexList,
					 list<ITL_trianglepatch>* lastTriangleList,
					 list<VECTOR3>* curVertexList,
					 list<ITL_trianglepatch>* curTriangleList,
					 int lastLevel,
					 int nDimension,
					 int nVerticesInATriangle,
					 int nNewTrianglesPerTriangle )
	{
		list<ITL_trianglepatch>::iterator iter;
		list<VECTOR3>::iterator iter2;
		int oldAID, oldBID, oldCID, a, b, c;
		VECTOR3 oldA, oldB, oldC;
		VECTOR3 ABMid, BCMid, CAMid;
		// a is ID of BCMid
		// b is ID of CAMid
		// c is ID of ABMid
		//int i = 0;
		int lastTriangleID = 0;
		int curVertexID = 0;
		double scaling;

		// Number of vertices and triangles in the last level
		int nVLast = (int)lastVertexList->size();
		int nTrLast = (int)lastTriangleList->size();
		#ifdef DEBUG_MODE
		printf( "Number of vertices and triangles in the last level: %d, %d\n", nVLast, nTrLast );
		#endif

		// First copy all the vertices of last level to the next level
		// Update vertex id
		curVertexList->insert( curVertexList->end(), lastVertexList->begin(), lastVertexList->end() );
		curVertexID = nVLast;

		// Iterate through each existing triangle and
		// create new triangles from it
		for( iter = lastTriangleList->begin(); iter != lastTriangleList->end(); iter++ )
		{
			// Get the vertex IDs of the existing triangle
			oldAID = (*iter).getVertexID( 0 );
			oldBID = (*iter).getVertexID( 1 );
			oldCID = (*iter).getVertexID( 2 );

			// Get the actual vertices (needed for computing mid point)
			iter2 = lastVertexList->begin();
			for( int k=0; k<oldAID; k++ )
				iter2++;
			oldA = (*iter2);

			iter2 = lastVertexList->begin();
			for( int k=0; k<oldBID; k++ )
				iter2++;
			oldB = (*iter2);

			iter2 = lastVertexList->begin();
			for( int k=0; k<oldCID; k++ )
				iter2++;
			oldC = (*iter2);

			// Compute the mid points of each side
			elementWiseMean( &oldA, &oldB, &ABMid, 1 );
			elementWiseMean( &oldB, &oldC, &BCMid, 1 );
			elementWiseMean( &oldC, &oldA, &CAMid, 1 );

			// Project the mid points on the spherical plane
			scaling = 1.0 / sqrt( ABMid.x()*ABMid.x() + ABMid.y()*ABMid.y() + ABMid.z()*ABMid.z() );
			ABMid.Set( ABMid.x()*scaling, ABMid.y()*scaling, ABMid.z()*scaling );
			scaling = 1.0 / sqrt( BCMid.x()*BCMid.x() + BCMid.y()*BCMid.y() + BCMid.z()*BCMid.z() );
			BCMid.Set( BCMid.x()*scaling, BCMid.y()*scaling, BCMid.z()*scaling );
			scaling = 1.0 / sqrt( CAMid.x()*CAMid.x() + CAMid.y()*CAMid.y() + CAMid.z()*CAMid.z() );
			CAMid.Set( CAMid.x()*scaling, CAMid.y()*scaling, CAMid.z()*scaling );

			// Add new vertices to the vertex list
			curVertexList->push_back( ABMid );
			curVertexList->push_back( BCMid );
			curVertexList->push_back( CAMid );
			int c = curVertexID;
			int a = curVertexID+1;
			int b = curVertexID+2;

			// Child Triangle 1
			ITL_trianglepatch nextPatch1( oldAID, c, b, lastLevel+1, lastTriangleID );
			curTriangleList->push_back( nextPatch1 );

			// Child Triangle 2
			ITL_trianglepatch nextPatch2( oldBID, c, a, lastLevel+1, lastTriangleID );
			curTriangleList->push_back( nextPatch2 );

			// Child Triangle 3
			ITL_trianglepatch nextPatch3( oldCID, a, b, lastLevel+1, lastTriangleID );
			curTriangleList->push_back( nextPatch3 );

			// Child Triangle 4
			ITL_trianglepatch nextPatch4( a, b, c, lastLevel+1, lastTriangleID );
			curTriangleList->push_back( nextPatch4 );

			// Update vertex and triangle IDs
			curVertexID = curVertexID + 3;
			lastTriangleID ++;

		}// end for

	}// end function

	// Average two arrays of vectors element-wise.
	// Each element of averages is the arithmetic mean of the
	// corresponding elements of the input arrays.
	void
	elementWiseMean( VECTOR3* A, VECTOR3* B, VECTOR3* avg, int nel )
	{
		assert( avg != NULL );
		for( int i=0; i<nel; i++ )
		{
			avg[i][0] = ( A[i].x() + B[i].x() ) / 2.0f;
			avg[i][1] = ( A[i].y() + B[i].y() ) / 2.0f;
			avg[i][2] = ( A[i].z() + B[i].z() ) / 2.0f;
		}

	}// end function

	// Computes the element-wise logarithms of the given values.
	// bases are the bases for which the logarithms are taken.
	double
	logarithm( double value, double base )
	{
		// This is a logarithmic identity.
		double l = log( value ) / log( base );
		return l;

	}// end function

	void
	computeTable( int nBin )
	{
		iNrOfThetas = nBin*2;							/**< Histogram parameter. Number of thetas.  */
		iNrOfPhis = nBin;								/**< Histogram parameter. Number of phis.  */
		fNrOfThetas = (float)iNrOfThetas;				/**< Histogram parameter. Number of phis.  */
		fNrOfPhis = (float)iNrOfPhis;					/**< Histogram parameter. Number of phis.  */
		piAngleMap = new int[iNrOfThetas*iNrOfPhis];	/**< Histogram parameter. Mapping from vector to a region in sperical coordinates.  */

		for( int t = 0; t < iNrOfThetas; t++ )
			for( int p = 0; p < iNrOfPhis; p++ )
			{
				float fTheta =	M_PI * 2.0f * float(t) / fNrOfThetas;
				float fPhi =	M_PI * float(p) / fNrOfPhis;
				int iBin = getBinByAngle( fTheta, fPhi, nBin );
				//cout << iBin << endl;
				if( iBin >= 0 )
					piAngleMap[t*iNrOfPhis + p] = iBin;
			}

		// set the flag to true
		bIsAngleMapInitialized = true;

	}// end function

	// Read the lookup table provided in the header file
	void
	computeTable1( int nBin, int nDivision )
	{
		list<ITL_trianglepatch>::iterator iter;
		list<VECTOR3>::iterator iter2;
		int AID, BID, CID;
		VECTOR3 A, B, C;

		// Compute area of triangle patch
		double triangleArea = (4*pi*1*1) / (double)nBin;

		// Compute allowable maximum length of table quad
		double quadLen = sqrt( triangleArea );

		// Compute minimum resolution of the table
		iNrOfThetas = (int)((2*pi)/(double)quadLen);
		iNrOfPhis = (int)(pi/(double)quadLen);
		fNrOfThetas = (float)iNrOfThetas;
		fNrOfPhis = (float)iNrOfPhis;

		//#ifdef DEBUG_MODE
		fprintf( stderr, "Triangular patch area: %g\n", triangleArea );
		fprintf( stderr, "Maximum allowable length of quad patch: %g\n", quadLen );
		fprintf( stderr, "Resolution of lookup table: %d %d\n", iNrOfThetas, iNrOfPhis );
		//#endif

		// Allocate memory for lookup table entries
		piAngleMap = new int[iNrOfThetas*iNrOfPhis];

		int unResolvedCount = 0;
		int index = 0;
		for( int t = 0; t < iNrOfThetas; t++ )
		{
			for( int p = 0; p < iNrOfPhis; p++ )
			{
				float fTheta =	M_PI * 2.0f * float(t) / fNrOfThetas;
				float fPhi =	M_PI * float(p) / fNrOfPhis;
				if( fTheta > 2*pi )	fTheta = 2*pi;
				if( fPhi > pi )	fPhi = pi;

				VECTOR3 v = spherical2Cartesian( 1, fTheta, fPhi );
				//cout << v[0] << " " << v[1] << " " << v[2] << endl;
				//cout << sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ) << endl;
				//v.Normalize();

				//int iBin = getBinNumber3D( v, &vertexList[nDivision-1], &triangleList[nDivision-1] );
				assert(  vertexList2[nDivision-1] != NULL );
				int iBin = getBinNumber3D( v, vertexList2[nDivision-1], &triangleList[nDivision-1] );

				if( iBin >= 0 )
					piAngleMap[t*iNrOfPhis + p] = iBin;
				else
				{
					unResolvedCount ++;
					printf( "Negative bin id: %d, theta: %g phi: %g\n", iBin, fTheta, fPhi );
				}

				//cout << index << endl;
				index ++;
			}// end inner for
		}// end outer for

		//printf( "Number of unresolved entries in the table: %d\n", unResolvedCount );

		// set the flag to true
		bIsAngleMapInitialized = true;
	}

	int
	//getBinNumber3D( VECTOR3 v, list<VECTOR3>* vertexList, list<ITL_trianglepatch>* triangleList )
	getBinNumber3D( VECTOR3 v, VECTOR3* vertexList, list<ITL_trianglepatch>* triangleList )
	{
		list<ITL_trianglepatch>::iterator iter;
		list<VECTOR3>::iterator iter2;
		int AID, BID, CID;
		VECTOR3 A, B, C, intersection;
		int binID;
		bool isIntersecting;

		//cout << "1"	<< endl;
		binID = 0;
		for( iter = triangleList->begin(); iter != triangleList->end(); iter++ )
		{
			// Get the vertex IDs of the existing triangle
			AID = (*iter).getVertexID( 0 );
			BID = (*iter).getVertexID( 1 );
			CID = (*iter).getVertexID( 2 );
			//printf( "%d %d %d\n", AID, BID, CID );

			//cout << "2"	<< endl;
			// Get the actual vertices (needed for computing mid point)
			//iter2 = vertexList->begin();
			//for( int k=0; k<AID; k++ )
			//	iter2++;
			//A = (*iter2);
			A = vertexList[AID];

			//cout << "3"	<< endl;
			//iter2 = vertexList->begin();
			//for( int k=0; k<BID; k++ )
			//	iter2++;
			//B = (*iter2);
			B = vertexList[BID];

			//cout << "4"	<< endl;
			//iter2 = vertexList->begin();
			//for( int k=0; k<CID; k++ )
			//	iter2++;
			//C = (*iter2);
			C = vertexList[CID];

			//cout << "5"	<< endl;

			isIntersecting = rayTriangleIntersection2( v, A, B, C, &intersection );
			if( isIntersecting )
			{
				return binID;
			}

			binID ++;
		}// end for

		//#ifdef DEBUG_MODE
		printf( "Vector falls in no bin ... \n" );
		//#endif
		return -1;

	}// end function

	int
	getBinNumber3DViaTable( VECTOR3 v, int nBin )
	{
		// convert the vector from Cartesian coodinates to sphercial coordinates;
		// the magnitude is ignored.
		float mytheta = getAngle(v.x(), v.y());//0~2pi
		float myphi = getAngle2(sqrt(v.x()*v.x()+v.y()*v.y()), v.z());//0~pi

		//if( !theta || !phi )
		//{
		//	printf("read in the regions first\n");
		//	return -1;
		//}

		if( bIsAngleMapInitialized == false )
			computeTable( nBin );

		// map the angle to the bin number
		int iTheta	=	min( iNrOfThetas - 1,	int( fNrOfThetas * mytheta / (M_PI * 2.0)));
		int iPhi	=	min( iNrOfPhis - 1, int( fNrOfPhis * myphi	/ M_PI));

		return piAngleMap[ iTheta * iNrOfPhis + iPhi];
	}


	// convert the angles in the spherical coordinates (mytheta, myphi) to the patch index
	// according to the lookup table composed of theta and phi.
	// The #entries in the lookup table is specified by the parameter 'binnum'.
	int
	getBinByAngle( float mytheta, float myphi, int binnum )
	{
		for(int i=0; i<binnum;i++)
		{
			if( (mytheta>=theta[i*2+0]) && (mytheta<=theta[i*2+1])&&
				(myphi>=phi[i*2+0]) && (myphi<=phi[i*2+1])
				)
			{
				return i;
			}
		}
		for(int i=0; i<binnum;i++)
		{
			if( ((mytheta+2*pi)>=theta[i*2+0]) && ((mytheta+2*pi)<=theta[i*2+1])&&
				(myphi>=phi[i*2+0]) && (myphi<=phi[i*2+1])
				)
			{
				return i;
			}
		}
		for(int i=0; i<binnum;i++)
		{
			if( ((mytheta)>=theta[i*2+0]) && ((mytheta)<=theta[i*2+1])&&
				((myphi+2*pi)>=phi[i*2+0]) && ((myphi+2*pi)<=phi[i*2+1] )
				)
			{
				return i;
			}
		}
		for(int i=0; i<binnum;i++)
		{
			if( ((mytheta+2*pi)>=theta[i*2+0]) && ((mytheta+2*pi)<=theta[i*2+1])&&
				((myphi+2*pi)>=phi[i*2+0]) && ((myphi+2*pi)<=phi[i*2+1])
				)
			{
				return i;
			}
		}
		return -1;
	}

	/*
	 * Assumption: ray originating at origin
	 */
	bool
	rayTriangleIntersection2( VECTOR3 d, VECTOR3 p0, VECTOR3 p1, VECTOR3 p2, VECTOR3* intersection )
	{
		double epsilon = 0.00001;

		//d.Normalize();

		// Normal to the triangle plane
		VECTOR3 e0( p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] );
		VECTOR3 e1( p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] );
		VECTOR3 q = cross( d, e1 );
		double a = dot( e0, q );

		if( a > (-epsilon) && a < epsilon )
		{
			(*intersection).Set( 0, 0, 0 );
			return false;
		}

		double f = 1.0 / a;
		VECTOR3 s( -p0[0], -p0[1], -p0[2] );
		double u = f * dot( s, q );

		if( u < 0 )
		{
			// Intersection is outside triangle
			(*intersection).Set( 0, 0, 0 );
			return false;
		}

		VECTOR3 r = cross( s, e0 );
		double v = f*dot( d, r );

		if( v<0 || (u+v)>1.0 )
		{
			// Intersection is outside triangle
			(*intersection).Set( 0, 0, 0 );
			return false;
		}

		double t = f * dot( e1, r );
		(*intersection).Set( t*d[0], t*d[1], t*d[2] );

		return true;


	}// end function

	void
	cartesian2Spherical( VECTOR3 v, double* mytheta, double* myphi )
	{
		(*mytheta) = getAngle(v.x(), v.y());//0~2pi
		(*myphi) =  getAngle2(sqrt(v.x()*v.x()+v.y()*v.y()), v.z());//0~pi
	}

	VECTOR3
	spherical2Cartesian( double r, double theta, double phi )
	{
		VECTOR3 ret;

		ret.Set( r * cos( theta ) * sin( phi ),
				 r * sin( theta ) * sin( phi ),
				 r * cos( phi ) );

		return ret;
	}

	double
	dot( VECTOR3 a, VECTOR3 b )
	{
		return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
	}

	VECTOR3
	cross( VECTOR3 a, VECTOR3 b )
	{
		return VECTOR3( a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] );
	}

	// compute the theta
	float
	getAngle(float x, float y)
	{

		if((x==0)&&(y==0))
			return 0;
		else
		{
			return (pi+(atan2(y,x)));
		}
	}

	// compute phi
	float
	getAngle2(float x, float y)
	{

		if((x==0)&&(y==0))
			return 0;
		else
		{
			return fabs(pi/2.0-(atan2(y,x)));
		}
	}

};

#endif
/* ITL_HISTOGRAMMAPPER_H_ */
