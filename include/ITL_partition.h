/**
 *  Class for space partitioning algorithms.
 *  A class that implements various space-partitioning algorithms.
 *  Created on: July 17, 2011.
 *  @author Abon
 *  @author Teng-Yok
 */
#ifndef ITL_PARTITION_H
#define ITL_PARTITION_H

#include "ITL_header.h"
#include "ITL_field_regular.h"
#include "ITL_globalentropy.h"

/*
float
distFeature( feature_t *F1, feature_t *F2 )
{
	//int dX = F1->X - F2->X, dY = F1->Y - F2->Y, dZ = F1->Z - F2->Z;
	//sqrt(dX*dX + dY*dY + dZ*dZ); 
	return abs( F1->binID - F2->binID );
}// End function
float (*funcPtr)( feature_t*, feature_t* ) = &distFeature;
*/

template <class T>
class ITL_partition
{
	int nDim;
	int mappingType;
	int nTrial;
	float percentArray[30];

public:
	/**
	 * Default constructor
	 */	
	ITL_partition( int nd, int mappingtype )
	{
		nDim = nd;
		mappingType = mappingtype;

		// Set the number of trials
		// and trial split positions
		nTrial = 13;
		percentArray[0] = 0.2f;
		percentArray[1] = 0.25f;
		percentArray[2] = 0.3f;
		percentArray[3] = 0.35f;
		percentArray[4] = 0.4f;
		percentArray[5] = 0.45f;
		percentArray[6] = 0.5f;
		percentArray[7] = 0.55f;
		percentArray[8] = 0.6f;
		percentArray[9] = 0.65f;
		percentArray[10] = 0.7f;
		percentArray[11] = 0.75f;
		percentArray[12] = 0.8f;

	}// End constructor

	/**
	 * Default destructor
	 */	
	~ITL_partition()
	{
	}

	/**
	 *
	*/
	void partition_Octree( float *parentLow, float *parentHigh,
						   ITL_field_regular<T> *parentField,
						   int nchild, int nBin,
						   float **childLow, float **childHigh,
						   float **childFreqList, float *entropyArray )
	{
		float mid[4];

		for( int i=0; i<nDim; i++ )
			mid[i] = floor( ( parentLow[i] + parentHigh[i] ) / 2.0f );

		// 1st child
		childLow[0][0] = parentLow[0]; childLow[0][1] = parentLow[1]; childLow[0][2] = parentLow[2];
		childHigh[0][0] = mid[0]; childHigh[0][1] = mid[1]; childHigh[0][2] = mid[2]; 
		// 2nd child
		childLow[1][0] = mid[0]; childLow[1][1] = parentLow[1]; childLow[1][2] = parentLow[2];
		childHigh[1][0] = parentHigh[0]; childHigh[1][1] = mid[1]; childHigh[1][2] = mid[2]; 
		// 3rd child
		childLow[2][0] = mid[0]; childLow[2][1] = mid[1]; childLow[2][2] = parentLow[2];
		childHigh[2][0] = parentHigh[0]; childHigh[2][1] = parentHigh[1]; childHigh[2][2] = mid[2]; 
		// 4th child
		childLow[3][0] = parentLow[0]; childLow[3][1] = mid[1]; childLow[3][2] = parentLow[2];
		childHigh[3][0] = mid[0]; childHigh[3][1] = parentHigh[1]; childHigh[3][2] = mid[2]; 
		// 5th child
		childLow[4][0] = parentLow[0]; childLow[4][1] = parentLow[1]; childLow[4][2] = mid[2];
		childHigh[4][0] = mid[0]; childHigh[4][1] = mid[1]; childHigh[4][2] = parentHigh[2]; 
		// 6th child
		childLow[5][0] = mid[0]; childLow[5][1] = parentLow[1]; childLow[5][2] = mid[2];
		childHigh[5][0] = parentHigh[0]; childHigh[5][1] = mid[1]; childHigh[5][2] = parentHigh[2]; 
		// 7th child
		childLow[6][0] = mid[0]; childLow[6][1] = mid[1]; childLow[6][2] = mid[2];
		childHigh[6][0] = parentHigh[0]; childHigh[6][1] = parentHigh[1]; childHigh[6][2] = parentHigh[2]; 
		// 8th child
		childLow[7][0] = parentLow[0]; childLow[7][1] = mid[1]; childLow[7][2] = mid[2];
		childHigh[7][0] = mid[0]; childHigh[7][1] = parentHigh[1]; childHigh[7][2] = parentHigh[2]; 	

		// Compute global entropy of each child 
		computeEntropyChildren( parentField, nchild, nBin,
								childLow, childHigh,
								childFreqList, entropyArray );	

	}// end function

	/**
	 *
	*/
	void partition_Quadtree( float *parentLow, float *parentHigh,
							 ITL_field_regular<T> *parentField,
							 int nchild, int nBin,
							 float **childLow, float **childHigh,
							 float **childFreqList, float *entropyArray )
	{
		float mid[4];
	
		for( int i=0; i<2; i++ )
			mid[i] = floor( ( parentLow[i] + parentHigh[i] ) / 2.0f );

		// 1st child
		childLow[0][0] = parentLow[0]; childLow[0][1] = parentLow[1]; childLow[0][2] = parentLow[2];
		childHigh[0][0] = mid[0]; childHigh[0][1] = mid[1]; childHigh[0][2] = parentHigh[2]; 
		// 2nd child
		childLow[1][0] = mid[0]; childLow[1][1] = parentLow[1]; childLow[1][2] = parentLow[2];
		childHigh[1][0] = parentHigh[0]; childHigh[1][1] = mid[1]; childHigh[1][2] = parentHigh[2]; 
		// 3rd child
		childLow[2][0] = mid[0]; childLow[2][1] = mid[1]; childLow[2][2] = parentLow[2];
		childHigh[2][0] = parentHigh[0]; childHigh[2][1] = parentHigh[1]; childHigh[2][2] = parentHigh[2]; 
		// 4th child
		childLow[3][0] = parentLow[0]; childLow[3][1] = mid[1]; childLow[3][2] = parentLow[2];
		childHigh[3][0] = mid[0]; childHigh[3][1] = parentHigh[1]; childHigh[3][2] = parentHigh[2]; 

		// Compute global entropy of each child 
		computeEntropyChildren( parentField, nchild, nBin,
								childLow, childHigh,
								childFreqList, entropyArray );	

	}// end function


	/**
	 *
	*/
	void partition_Kdtree_SplitMiddle( float *parentLow, float *parentHigh,
									   ITL_field_regular<T> *parentField,
									   int level, int nchild, int nBin,
									   float **childLow, float **childHigh,
									   float **childFreqList, float *entropyArray )
	{
		int selectedDim = level % nDim;
		float mid = floor( ( parentLow[selectedDim] + parentHigh[selectedDim] ) / 2.0f );

		for( int i = 0; i<nDim; i++ )
		{
			// Create partitions
			if( i == selectedDim )
			{
				childLow[0][i] = parentLow[i]; childHigh[0][i] = mid;
				childLow[1][i] = mid; childHigh[1][i] = parentHigh[i];
			}
			else
			{
				childLow[0][i] = parentLow[i]; childHigh[0][i] = parentHigh[i];
				childLow[1][i] = parentLow[i]; childHigh[1][i] = parentHigh[i];
			}
		}// end for

		// Compute global entropy of each child 
		computeEntropyChildren( parentField, nchild, nBin,
								childLow, childHigh,
								childFreqList, entropyArray );	
			
	}// end function

	/**
	 *
	*/
	void partition_Kdtree_SplitMiddle2( float *parentLow, float *parentHigh,
										ITL_field_regular<int> *parentBinField,
										int level, int nchild, int nBin,
										float **childLow, float **childHigh,
										float **childFreqList, float *entropyArray )
	{
		int selectedDim = level % 2;
		float mid = floor( ( parentLow[selectedDim] + parentHigh[selectedDim] ) / 2.0f );

		for( int i = 0; i<nDim; i++ )
		{
			// Create partitions
			if( i == selectedDim )
			{
				childLow[0][i] = parentLow[i]; childHigh[0][i] = mid;
				childLow[1][i] = mid; childHigh[1][i] = parentHigh[i];
			}
			else
			{
				childLow[0][i] = parentLow[i]; childHigh[0][i] = parentHigh[i];
				childLow[1][i] = parentLow[i]; childHigh[1][i] = parentHigh[i];
			}
		}// end for

		// Compute global entropy of each child 
		computeEntropyChildren2( parentBinField, nchild, nBin,
								 childLow, childHigh,
								 childFreqList, entropyArray );	
		
	
	}// end function

	/**
	 *
	 */
	void partition_MaxEntropy( float *parentLow, float *parentHigh,
							   ITL_field_regular<T> *parentField,
							   int nchild, int nBin,
							   float **childLow, float **childHigh,
							   float **childFreqList, float *entropyArray,
							   FILE* dumpFile = NULL )
	{	
		// Initialize local variables		
		float mid = 0;		
		#if defined( _WIN32 ) || defined( _WIN64 )
			float* entropyArrayCopy = new float[nchild];
		#endif
		#if defined( _LINUX )
			float entropyArrayCopy[nchild];
		#endif
		
		float scoreCopy = -100000;
		for( int i=0; i<nchild; i++ )
			entropyArrayCopy[i] = -100000;
		float **childLowCopy = new float*[2];
		float **childHighCopy = new float*[2];
		for( int i=0; i<nchild; i++ )
		{
			childLowCopy[i] = new float[nDim];
			childHighCopy[i] = new float[nDim];
		}

		fprintf( stderr, "Parent: (%f %f %f) -- (%f %f %f)\n",
			parentLow[0], parentLow[1], parentLow[2],
			parentHigh[0], parentHigh[1], parentHigh[2] );
		if( dumpFile != NULL )
			fprintf( dumpFile, "%f, %f, %f, %f, %f, %f\n",
				parentLow[0], parentLow[1], parentLow[2],
				parentHigh[0], parentHigh[1], parentHigh[2] );
				
		for( int i = 0; i<nDim-1; i++ )		
		//for( int i = 0; i<nDim; i++ )
		{
			for( int j = 0; j<nTrial; j++ )
			{
				mid = floor( parentLow[i] + percentArray[j]*( parentHigh[i]-parentLow[i] ) );

				// Create partitions
				for( int k=0; k<nDim; k++ )
				{
					if( k == i )
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = mid;
						childLow[1][k] = mid; childHigh[1][k] = parentHigh[k];
					}
					else
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = parentHigh[k];
						childLow[1][k] = parentLow[k]; childHigh[1][k] = parentHigh[k];

					}
				}

				printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n",
					childLow[0][0], childLow[0][1], childLow[0][2],
					childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n",
					childLow[1][0], childLow[1][1], childLow[1][2],
					childHigh[1][0], childHigh[1][1], childHigh[1][2] );

				// Compute global entropy of each child 
				computeEntropyChildren( parentField, nchild, nBin,
										childLow, childHigh,
										childFreqList, entropyArray );	
				
				fprintf( stderr, "%f, %f\n", entropyArray[0], entropyArray[1] );
				if( dumpFile != NULL )
					fprintf( dumpFile, "%f, %f\n", entropyArray[0], entropyArray[1] );

				// Remember this partition, if it maximizes entropy of the two children				
				//if( max( entropyArray[0], entropyArray[1] ) > max( entropyArrayCopy[0], entropyArrayCopy[1] ) )
				float score = max( entropyArray[0], entropyArray[1] );
				//float score = ( entropyArray[0]*percntArray[j]+entropyArray[1]*(1-percntArray[j]) );
				//float score = fabs( entropyArray[0]*percntArray[j] - entropyArray[1]*(1-percntArray[j]) );
				//float t1 = ( entropyArray[0] / ( (float)log( (float)nBin ) / log(2.0f ) ) ) * percntArray[j];
 				//float t2 = ( entropyArray[1] / ( (float)log( (float)nBin ) / log(2.0f ) ) ) * (1-percntArray[j]);
				//float score = (t1 + t2)/2.0f + fabs( t1-t2 );
				//float score = ( entropyArray[0]+entropyArray[1] ) * fabs( entropyArray[0]-entropyArray[1]  );
				if( score > scoreCopy )
				{
					//cout << "in" << endl;
					for( int iC=0; iC<nchild; iC++ )
					{
						entropyArrayCopy[iC] = entropyArray[iC];						
						memcpy( childLowCopy[iC], childLow[iC], sizeof(float)*nDim );
						memcpy( childHighCopy[iC], childHigh[iC], sizeof(float)*nDim );						
					}
					scoreCopy = score;
				}

			}// end for j: trials
		}// end for i: dimensions

		// Return the stored partition which maximizes entropy
		for( int iC=0; iC<nchild; iC++ )
		{
			entropyArray[iC] = entropyArrayCopy[iC];	
			memcpy( childLow[iC], childLowCopy[iC], sizeof(float)*nDim );
			memcpy( childHigh[iC], childHighCopy[iC], sizeof(float)*nDim );						
		}
		
		for( int i=0; i<nchild; i++ )
		{
			delete [] childLowCopy[i];
			delete [] childHighCopy[i];
		}
		#if defined( _WIN32 ) || defined( _WIN64 )
			delete [] entropyArrayCopy;
		#endif

	
	}// end function


	/**
	 *
	 */
	void partition_MaxEntropy2( float *parentLow, float *parentHigh,
								ITL_field_regular<int> *parentBinField,
								int nchild, int nBin,
								float **childLow, float **childHigh,
								float **childFreqList, float *entropyArray )
	{	
		// Initialize local variables		
		float mid = 0;	
		float entropyArrayCopy[nchild];
		for( int i=0; i<nchild; i++ )
			entropyArrayCopy[i] = -100000;
		float **childLowCopy = new float*[2];
		float **childHighCopy = new float*[2];
		for( int i=0; i<nchild; i++ )
		{
			childLowCopy[i] = new float[nDim];
			childHighCopy[i] = new float[nDim];
		}

		//for( int i = 0; i<nDim; i++ )		
		for( int i = 0; i<nDim-1; i++ )
		{
			for( int j = 0; j<nTrial; j++ )
			{
				mid = floor( parentLow[i] + percntArray[j]*( parentHigh[i]-parentLow[i] ) );

				// Create partitions
				for( int k=0; k<nDim; k++ )
				{
					if( k == i )
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = mid;
						childLow[1][k] = mid; childHigh[1][k] = parentHigh[k];
					}
					else
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = parentHigh[k];
						childLow[1][k] = parentLow[k]; childHigh[1][k] = parentHigh[k];

					}
				}

				/*
				printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n",
					childLow[0][0], childLow[0][1], childLow[0][2],
					childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n",
					childLow[1][0], childLow[1][1], childLow[1][2],
					childHigh[1][0], childHigh[1][1], childHigh[1][2] );
				*/

				// Compute global entropy of each child 
				computeEntropyChildren2( parentBinField, nchild, nBin,
										 childLow, childHigh,
										 childFreqList, entropyArray );	
				//cout << entropyArray[0] << " " << entropyArray[1] << endl;
				//cout << entropyArrayCopy[0] << " " << entropyArrayCopy[1] << endl;

				// Remember this partition, if it maximizes entropy of the two children				
				//if( max( entropyArray[0], entropyArray[1] ) >
				//    max( entropyArrayCopy[0], entropyArrayCopy[1] ) )
				if( ( entropyArray[0]+entropyArray[1] ) >
					( entropyArrayCopy[0]+entropyArrayCopy[1] ) )
				{
					//cout << "in" << endl;
					for( int iC=0; iC<nchild; iC++ )
					{
						entropyArrayCopy[iC] = entropyArray[iC];						
						memcpy( childLowCopy[iC], childLow[iC], sizeof(float)*nDim );
						memcpy( childHighCopy[iC], childHigh[iC], sizeof(float)*nDim );						
					}
				}


			}// end for j: trials
		}// end for i: dimensions

		// Return the stored partition which maximizes entropy
		for( int iC=0; iC<nchild; iC++ )
		{
			entropyArray[iC] = entropyArrayCopy[iC];	
			memcpy( childLow[iC], childLowCopy[iC], sizeof(float)*nDim );
			memcpy( childHigh[iC], childHighCopy[iC], sizeof(float)*nDim );						
		}
		
		for( int i=0; i<nchild; i++ )
		{
			delete [] childLowCopy[i];
			delete [] childHighCopy[i];
		}

	}// end function

	/**
	 * EMD based space-partitioning function.
	 * This function partitions a block in such a way
	 * that each child has distributions close to that
	 * of the maximum entropy distribution. This 
	 * particular version passes down the actual data
	 * and hence consumes more memory.
	 */
	int
	partition_EMD( float *parentLow, float *parentHigh,
				   ITL_field_regular<T> *parentField,
				   int nchild, int nBin,
				   float **childLow, float **childHigh,
				   float **childFreqList, float *metricArray,
				   FILE* dumpFile = NULL )
	{	
		// Initialize local variables		
		float mid = 0;

		// If the field has very low z-thickness
		// then only partition along x and y
		int upperDim;
		if( mappingType == 2 )
			upperDim = 2;
		else if( mappingType == 3 )
			upperDim = 3;

		#if defined( _WIN32 ) || defined( _WIN64 )
			float* tempMetricArrayLeft = new float[nTrial];							
			float* tempMetricArrayRight = new float[nTrial];							
			float* diffMetricArrayLeft = new float[nTrial];							
			float* diffMetricArrayRight = new float[nTrial];							
			int* diffSignArrayLeft = new int[nTrial];							
			int* diffSignArrayRight = new int[nTrial];	
		#endif
		#if defined( _LINUX )
			float tempMetricArrayLeft[nTrial];
			float tempMetricArrayRight[nTrial];
			float diffMetricArrayLeft[nTrial];
			float diffMetricArrayRight[nTrial];
			int diffSignArrayLeft;
			int diffSignArrayRight;
		#endif
		int nInfPointLeft, nInfPointRight;
		float metricArrayCopy[2];
		float scoreCopyArray[] = { 100000, 100000 } ;
		float **childLowCopy = new float*[2];
		float **childHighCopy = new float*[2];
		for( int i=0; i<nchild; i++ )
		{
			childLowCopy[i] = new float[nDim];
			childHighCopy[i] = new float[nDim];
		}
		int stopPartition = 0;

		//for( int i=0; i<nchild; i++ )
		//	metricArrayCopy[i] = -100000;
		//float **childLowCopy = new float*[2];
		//float **childHighCopy = new float*[2];
		//for( int i=0; i<nchild; i++ )
		//{
		//	childLowCopy[i] = new float[nDim];
		//	childHighCopy[i] = new float[nDim];
		//}

		fprintf( stderr, "Parent: (%f %f %f) -- (%f %f %f)\n",
			parentLow[0], parentLow[1], parentLow[2],
			parentHigh[0], parentHigh[1], parentHigh[2] );
		if( dumpFile != NULL )
			fprintf( dumpFile, "%f, %f, %f, %f, %f, %f\n",
				parentLow[0], parentLow[1], parentLow[2],
				parentHigh[0], parentHigh[1], parentHigh[2] );
				
		// Iterate over each dimension
		// and for each dimension, iterate
		// over all trial split positions
		int trialIndex = 0;
		float tempArray[2];
		for( int i = 0; i<upperDim; i++ )		
		{
			trialIndex = 0;
			for( int j = 0; j<nTrial; j++ )
			{
				mid = floor( parentLow[i] + percntArray[j]*
							( parentHigh[i]-parentLow[i] ) );

				// Create partitions
				for( int k=0; k<nDim; k++ )
				{
					if( k == i )
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = mid;
						childLow[1][k] = mid; childHigh[1][k] = parentHigh[k];
					}
					else
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = parentHigh[k];
						childLow[1][k] = parentLow[k]; childHigh[1][k] = parentHigh[k];
					}
				}

				printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n",
					childLow[0][0], childLow[0][1], childLow[0][2],
					childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n",
					childLow[1][0], childLow[1][1], childLow[1][2],
					childHigh[1][0], childHigh[1][1], childHigh[1][2] );

				// Compute EMD based metric of each child 
				computeEMDMetricChildren( parentField, nchild, nBin,
										  childLow, childHigh,
										  childFreqList, tempArray );
				tempMetricArrayLeft[trialIndex] = tempArray[0];
				tempMetricArrayRight[trialIndex] = tempArray[1];
			
				// Write computed metric values
				fprintf( stderr, "%g, %g\n", tempMetricArrayLeft[trialIndex],
											 tempMetricArrayRight[trialIndex] );
				if( dumpFile != NULL )
					fprintf( dumpFile, "%g, %g\n", tempMetricArrayLeft[trialIndex],
												   tempMetricArrayRight[trialIndex] );

				trialIndex ++;

			}// end for j: trials

			// For this dimension,
			// find the split point where the metric
			// curve hits a low
			// Update split position based on first curve
			differentiate( tempMetricArrayLeft, diffMetricArrayLeft, nTrial );
			findMinPoints( diffMetricArrayLeft, nTrial,
						   diffSignArrayLeft, &nInfPointLeft  );
			updateSplitPosition( &stopPartition, i,
								 nInfPointLeft,
								 tempMetricArrayLeft,
								 diffSignArrayLeft,
								 parentLow, parentHigh,
								 childLowCopy, childHighCopy,
								 scoreCopyArray );

			// Update split position based on second curve
			differentiate( tempMetricArrayRight, diffMetricArrayRight, nTrial );
			findMinPoints( diffMetricArrayRight, nTrial,
						   diffSignArrayRight, &nInfPointRight );
			updateSplitPosition( &stopPartition, i,
								 nInfPointRight,
								 tempMetricArrayRight,
								 diffSignArrayRight,
								 parentLow, parentHigh,
								 childLowCopy, childHighCopy,
								 scoreCopyArray );
	
		}// end for i: dimensions

		// Return the stored partition which 
		// minimized EMD based metric
		if( stopPartition == 0 )
		{
			for( int iC=0; iC<nchild; iC++ )
			{
				metricArray[iC] = scoreCopyArray[iC];	
				memcpy( childLow[iC], childLowCopy[iC], sizeof(float)*nDim );
				memcpy( childHigh[iC], childHighCopy[iC], sizeof(float)*nDim );						
			}
			printf( "Selected child 1: (%f %f %f) -- (%f %f %f)\n",
				childLow[0][0], childLow[0][1], childLow[0][2],
				childHigh[0][0], childHigh[0][1], childHigh[0][2] );
			printf( "Selected child 2: (%f %f %f) -- (%f %f %f)\n",
				childLow[1][0], childLow[1][1], childLow[1][2],
				childHigh[1][0], childHigh[1][1], childHigh[1][2] );
		}
		
		for( int i=0; i<nchild; i++ )
		{
			delete [] childLowCopy[i];
			delete [] childHighCopy[i];
		}

		return stopPartition;
	
	}// end function

	/**
	 * EMD based space-partitioning function.
	 * This function partitions a block in such a way
	 * that each child has distributions close to that
	 * of the maximum entropy distribution. This 
	 * particular version passes down the histogram binIDs
	 * instead of the actual data.
	 */
	int
	partition_EMD2( float *parentLow, float *parentHigh,
				    ITL_field_regular<int> *parentBinField,
					int nchild, int nBin,
					float **childLow, float **childHigh,
					float **childFreqList, float *metricArray,
					FILE* dumpFile = NULL )
	{	
		// Initialize local variables		
		float mid = 0;

		// If the field has very low z-thickness
		// then only partition along x and y
		int upperDim;
		if( mappingType == 2 )
			upperDim = 2;
		else if( mappingType == 3 )
			upperDim = 3;

		#if defined( _WIN32 ) || defined( _WIN64 )
			float* tempMetricArrayLeft = new float[nTrial];							
			float* tempMetricArrayRight = new float[nTrial];							
			float* diffMetricArrayLeft = new float[nTrial];							
			float* diffMetricArrayRight = new float[nTrial];							
			int* diffSignArrayLeft = new int[nTrial];							
			int* diffSignArrayRight = new int[nTrial];	
		#endif
		#if defined( _LINUX )
			float tempMetricArrayLeft[nTrial];
			float tempMetricArrayRight[nTrial];
			float diffMetricArrayLeft[nTrial];
			float diffMetricArrayRight[nTrial];
			int diffSignArrayLeft;
			int diffSignArrayRight;
		#endif
		int nInfPointLeft, nInfPointRight;
		float scoreCopyArray[] = { 100000, 100000 } ;
		float **childLowCopy = new float*[2];
		float **childHighCopy = new float*[2];
		for( int i=0; i<nchild; i++ )
		{
			childLowCopy[i] = new float[nDim];
			childHighCopy[i] = new float[nDim];
		}
		int stopPartition = 0;

		//for( int i=0; i<nchild; i++ )
		//	metricArrayCopy[i] = -100000;
		//float **childLowCopy = new float*[2];
		//float **childHighCopy = new float*[2];
		//for( int i=0; i<nchild; i++ )
		//{
		//	childLowCopy[i] = new float[nDim];
		//	childHighCopy[i] = new float[nDim];
		//}

		fprintf( stderr, "Parent: (%f %f %f) -- (%f %f %f)\n",
			parentLow[0], parentLow[1], parentLow[2],
			parentHigh[0], parentHigh[1], parentHigh[2] );
		if( dumpFile != NULL )
			fprintf( dumpFile, "%f, %f, %f, %f, %f, %f\n",
				parentLow[0], parentLow[1], parentLow[2],
				parentHigh[0], parentHigh[1], parentHigh[2] );
				
		// Iterate over each dimension
		// and for each dimension, iterate
		// over all trial split positions
		int trialIndex = 0;
		float tempArray[2];
		// Outer for loop: each dimension
		for( int i = 0; i<upperDim; i++ )		
		{
			trialIndex = 0;
			// Inner for loop: each trial
			for( int j = 0; j<nTrial; j++ )
			{
				mid = floor( parentLow[i] + percentArray[j]*
							( parentHigh[i]-parentLow[i] ) );

				// Create partitions
				for( int k=0; k<nDim; k++ )
				{
					if( k == i )
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = mid;
						childLow[1][k] = mid; childHigh[1][k] = parentHigh[k];
					}
					else
					{
						childLow[0][k] = parentLow[k]; childHigh[0][k] = parentHigh[k];
						childLow[1][k] = parentLow[k]; childHigh[1][k] = parentHigh[k];
					}
				}

				printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n",
					childLow[0][0], childLow[0][1], childLow[0][2],
					childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n",
					childLow[1][0], childLow[1][1], childLow[1][2],
					childHigh[1][0], childHigh[1][1], childHigh[1][2] );

				// Compute EMD based metric of each child 
				computeEMDMetricChildren2( parentBinField, nchild, nBin,
										  childLow, childHigh,
										  childFreqList, tempArray );
				tempMetricArrayLeft[trialIndex] = tempArray[0];
				tempMetricArrayRight[trialIndex] = tempArray[1];
			
				// Write computed metric values
				fprintf( stderr, "%g, %g\n", tempMetricArrayLeft[trialIndex],
											 tempMetricArrayRight[trialIndex] );
				if( dumpFile != NULL )
					fprintf( dumpFile, "%g, %g\n", tempMetricArrayLeft[trialIndex],
												   tempMetricArrayRight[trialIndex] );

				trialIndex ++;

			}// end for j: trials

			// For this dimension,
			// find the split point where the metric
			// curve hits a low
			// Update split position based on first curve
			differentiate( tempMetricArrayLeft, diffMetricArrayLeft, nTrial );
			findMinPoints( diffMetricArrayLeft, nTrial,
						   diffSignArrayLeft, &nInfPointLeft  );
			updateSplitPosition( &stopPartition, i,
								 nInfPointLeft,
								 tempMetricArrayLeft,
								 diffSignArrayLeft,
								 parentLow, parentHigh,
								 childLowCopy, childHighCopy,
								 scoreCopyArray );

			// Update split position based on second curve
			differentiate( tempMetricArrayRight, diffMetricArrayRight, nTrial );
			findMinPoints( diffMetricArrayRight, nTrial,
						   diffSignArrayRight, &nInfPointRight );
			updateSplitPosition( &stopPartition, i,
								 nInfPointRight,
								 tempMetricArrayRight,
								 diffSignArrayRight,
								 parentLow, parentHigh,
								 childLowCopy, childHighCopy,
								 scoreCopyArray );
	
		}// end for i: dimensions

		// Return the stored partition which 
		// minimized EMD based metric
		if( stopPartition == 0 )
		{
			for( int iC=0; iC<nchild; iC++ )
			{
				metricArray[iC] = scoreCopyArray[iC];	
				memcpy( childLow[iC], childLowCopy[iC], sizeof(float)*nDim );
				memcpy( childHigh[iC], childHighCopy[iC], sizeof(float)*nDim );						
			}
			printf( "Selected child 1: (%f %f %f) -- (%f %f %f)\n",
				childLow[0][0], childLow[0][1], childLow[0][2],
				childHigh[0][0], childHigh[0][1], childHigh[0][2] );
			printf( "Selected child 2: (%f %f %f) -- (%f %f %f)\n",
				childLow[1][0], childLow[1][1], childLow[1][2],
				childHigh[1][0], childHigh[1][1], childHigh[1][2] );
		}
		
		for( int i=0; i<nchild; i++ )
		{
			delete [] childLowCopy[i];
			delete [] childHighCopy[i];
		}

		return stopPartition;
	
	}// End function

	void	
	updateSplitPosition( int *stopPartition,
						 int selectDim, int nInfPoint,
						 float *tempMetricArray,
						 int *diffSignArray,
						 float *parentLow, float *parentHigh,
						 float **childLowCopy, float **childHighCopy,
						 float *scoreCopyArray )
	{
		fprintf( stderr,
				 "Number of minima on EMD curve: %d\n",
				 nInfPoint );

		// Set stopping criterion
		if( nInfPoint == 0 )	(*stopPartition) = 1;
		else					(*stopPartition) = 0;
		
		// Scan through the minima on the curve
		// to find the optimal split position
		float mid;
		for( int iP=0; iP<nInfPoint; iP++ )
		{
			fprintf( stderr,
					 "Position of next minima"
					 "on EMD curve: %d\n", diffSignArray[iP] );
		
			// Update split position 
			// if following condition is satisfied
			if( tempMetricArray[diffSignArray[iP]] < 
				min( scoreCopyArray[0], scoreCopyArray[1] )  )
			{
				scoreCopyArray[0] = tempMetricArray[diffSignArray[iP]];
				scoreCopyArray[1] = tempMetricArray[diffSignArray[iP]];
				
				mid = floor( parentLow[selectDim] +
							 percentArray[diffSignArray[iP]]*
							( parentHigh[selectDim]-parentLow[selectDim] ) );

				// Create the partition
				// from parent block information
				for( int k=0; k<nDim; k++ )
				{
					if( k == selectDim )
					{
						childLowCopy[0][k] = parentLow[k];
						childHighCopy[0][k] = mid;
						childLowCopy[1][k] = mid;
						childHighCopy[1][k] = parentHigh[k];
					}
					else
					{
						childLowCopy[0][k] = parentLow[k];
						childHighCopy[0][k] = parentHigh[k];
						childLowCopy[1][k] = parentLow[k];
						childHighCopy[1][k] = parentHigh[k];
					}
				}// End inner for
			}// End if
		}// End outer for			

	}// End function

	void computeEntropyChildren( ITL_field_regular<T> *parentField,
								 int nchild, int nBin,
								 float **childLow, float **childHigh,
								 float **childFreqList, float *entropyArray )
	{
		//int lowPad[3];
		//int highPad[3];
		float lowF[4];
		float highF[4];
		//int sizeNeighborhoodArray[3];
		lowF[0] = lowF[1] = lowF[2] = 0;
		//lowPad[0] = lowPad[1] = lowPad[2] = 0;
		//highPad[0] = highPad[1] = highPad[2] = 0;
		//sizeNeighborhoodArray[0] = sizeNeighborhoodArray[1] = sizeNeighborhoodArray[2] = 0;
		T* childFieldData = NULL;
		ITL_field_regular<T> *childField = NULL;
		ITL_globalentropy<T> *globalEntropyComputer = NULL;

		for( int i=0; i<nchild; i++ )
		{
			// Get data for the next child
			//printf( "%f -- %f -- %f\n", childLow[i][0], childLow[i][1], childLow[i][2] );
			//printf( "%f -- %f -- %f\n", childHigh[i][0], childHigh[i][1], childHigh[i][2] );
			childFieldData = parentField->getDataBetween( childLow[i], childHigh[i] );

			// Create sub-field
			highF[0] = childHigh[i][0] - childLow[i][0];
			highF[1] = childHigh[i][1] - childLow[i][1];
			highF[2] = childHigh[i][2] - childLow[i][2];
			childField = new ITL_field_regular<T>( childFieldData, nDim, lowF, highF );
										//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute global entropy of the sub-field
			globalEntropyComputer = new ITL_globalentropy<T>( childField );
			if( mappingType == 2 )
				globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			else if( mappingType == 3 )
				globalEntropyComputer->computeHistogramBinField( "vector", nBin );			
			globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );
			entropyArray[i] = globalEntropyComputer->getGlobalEntropy();
			globalEntropyComputer->getHistogramFrequencies( nBin, childFreqList[i] );
			
			/*
			printf( "%f, %f, %f, %f, %f, %f\n ",
			childLow[i][0], childHigh[i][0], childLow[i][1],
			childHigh[i][1], childLow[i][2], childHigh[i][2] );   						
			for( int k=0; k<nBin; k++ )
				printf( "%d, ", childFreqList[i][k] );
			printf( "\n" );
			*/
			
			delete globalEntropyComputer;
			delete childField;
			delete childFieldData;

		}// end for

		
	}// end function

	void computeEntropyChildren2( ITL_field_regular<int> *parentBinField,
								  int nchild, int nBin,
								  float **childLow, float **childHigh,
								  float **childFreqList, float *entropyArray )
	{
		float lowF[3];
		float highF[3];
		lowF[0] = lowF[1] = lowF[2] = 0;
		int* childFieldBinData = NULL;
		ITL_field_regular<T> *childField = NULL;
		ITL_field_regular<int> *childBinField = NULL;
		ITL_globalentropy<T> *globalEntropyComputer = NULL;

		for( int i=0; i<nchild; i++ )
		{
			// Get histogram bin data for the next child
			childFieldBinData = parentBinField->getDataBetween( childLow[i], childHigh[i] );

			// Create sub-field
			highF[0] = childHigh[i][0] - childLow[i][0];
			highF[1] = childHigh[i][1] - childLow[i][1];
			highF[2] = childHigh[i][2] - childLow[i][2];

			childField = new ITL_field_regular<T>( nDim, lowF, highF );
			//, lowPad, highPad, sizeNeighborhoodArray );
			childBinField = new ITL_field_regular<int>( childFieldBinData, nDim, lowF, highF );
			//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute global entropy of the sub-field
			globalEntropyComputer = new ITL_globalentropy<T>( NULL, childBinField );
			if( mappingType == 2 )
				globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			else if( mappingType == 3 )
				globalEntropyComputer->computeHistogramBinField( "vector", nBin );			
			globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );

			entropyArray[i] = globalEntropyComputer->getGlobalEntropy();

			delete globalEntropyComputer;
			delete childField;
			//delete childBinField;

		}// end for
		
	}// end function

	//static	
	void
	computeEMDMetricChildren( ITL_field_regular<T> *parentField,
	 					      int nchild, int nBin,
							  float **childLow, float **childHigh,
							  float **childFreqList, float *metricArray )
	{
		float lowF[4];
		float highF[4];
		lowF[0] = lowF[1] = lowF[2] = 0;
		T* childFieldData = NULL;
		ITL_field_regular<T> *childField = NULL;
		ITL_globalentropy<T> *globalEntropyComputer = NULL;
		//float *uniformWeight = new float[nBin];
		float uniformWeight[360];
		float uWeight = 1.0f / (float)nBin;
		for( int i=0; i<nBin; i++ )
			uniformWeight[i] = uWeight;

		for( int i=0; i<nchild; i++ )
		{
			// Get data for the next child
			//printf( "%f -- %f -- %f\n", childLow[i][0], childLow[i][1], childLow[i][2] );
			//printf( "%f -- %f -- %f\n", childHigh[i][0], childHigh[i][1], childHigh[i][2] );
			childFieldData = parentField->getDataBetween( childLow[i], childHigh[i] );

			// Create sub-field
			highF[0] = childHigh[i][0] - childLow[i][0];
			highF[1] = childHigh[i][1] - childLow[i][1];
			highF[2] = childHigh[i][2] - childLow[i][2];
			childField = new ITL_field_regular<T>( childFieldData, nDim, lowF, highF );
										//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute EMD-based metric of sub field
			globalEntropyComputer = new ITL_globalentropy<T>( childField );
			if( mappingType == 2 )
				globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			else if( mappingType == 3 )
				globalEntropyComputer->computeHistogramBinField( "vector", nBin );			
			globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );
			globalEntropyComputer->getHistogramFrequencies( nBin, childFreqList[i] );
			metricArray[i] = computeEMD1D( uniformWeight, childFreqList[i], nBin );
			
			/*
			printf( "%f, %f, %f, %f, %f, %f\n ",
			childLow[i][0], childHigh[i][0], childLow[i][1],
			childHigh[i][1], childLow[i][2], childHigh[i][2] );   						
			for( int k=0; k<nBin; k++ )
				printf( "%d, ", childFreqList[i][k] );
			printf( "\n" );
			*/
			
			delete globalEntropyComputer;
			delete childField;
			delete childFieldData;
			//delete [] uniformWeight;

		}// end for
		
	}// End function

	void
	computeEMDMetricChildren2( ITL_field_regular<int> *parentBinField,
	 					       int nchild, int nBin,
							   float **childLow, float **childHigh,
							   float **childFreqList, float *metricArray )
	{
		float lowF[4];
		float highF[4];
		lowF[0] = lowF[1] = lowF[2] = 0;
		//T* childFieldData = NULL;
		//ITL_field_regular<T> *childField = NULL;
		//ITL_field_regular<int> *childBinField = NULL;
		ITL_globalentropy<T> *globalEntropyComputer = NULL;
		//float *uniformWeight = new float[nBin];
		float uniformWeight[360];
		float uWeight = 1.0f / (float)nBin;
		for( int i=0; i<nBin; i++ )
			uniformWeight[i] = uWeight;

		int *childBinFieldData =NULL;
		for( int i=0; i<nchild; i++ )
		{
			// Get data for the next child
			//printf( "%f -- %f -- %f\n", childLow[i][0], childLow[i][1], childLow[i][2] );
			//printf( "%f -- %f -- %f\n", childHigh[i][0], childHigh[i][1], childHigh[i][2] );
			childBinFieldData = parentBinField->getDataBetween( childLow[i], childHigh[i] );

			// Create sub-field
			highF[0] = childHigh[i][0] - childLow[i][0];
			highF[1] = childHigh[i][1] - childLow[i][1];
			highF[2] = childHigh[i][2] - childLow[i][2];
			//childField = new ITL_field_regular<T>( NULL, nDim, lowF, highF );
			//childField->setBinfield( 
										//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute EMD-based metric of sub field
			//globalEntropyComputer = new ITL_globalentropy<T>( childField );
			//if( mappingType == 2 )
			//	globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			//else if( mappingType == 3 )
			//	globalEntropyComputer->computeHistogramBinField( "vector", nBin );			
			//globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );
			//globalEntropyComputer->getHistogramFrequencies( nBin, childFreqList[i] );
			int nPointChild = (highF[0]+1)*(highF[1]+1)*(highF[2]+1);
			computeHistogramFrequencies( childBinFieldData, nPointChild,
									     childFreqList[i], nBin );
			delete [] childBinFieldData;

			metricArray[i] = computeEMD1D( uniformWeight, childFreqList[i], nBin );
			
			/*
			printf( "%f, %f, %f, %f, %f, %f\n ",
			childLow[i][0], childHigh[i][0], childLow[i][1],
			childHigh[i][1], childLow[i][2], childHigh[i][2] );   						
			for( int k=0; k<nBin; k++ )
				printf( "%d, ", childFreqList[i][k] );
			printf( "\n" );
			*/
			
			//delete globalEntropyComputer;
			//delete childField;
			//delete childFieldData;
			//delete [] uniformWeight;

		}// end for
		
	}// End function

	void
	computeHistogramFrequencies( int *binIdList, int nPoint,
								 float *freqList, int nBin )
	{
		assert( freqList != NULL );
		for( int i=0; i<nBin; i++ )
			freqList[i] = 0;
		for( int i=0; i<nPoint; i++ )
			freqList[ binIdList[i] ] ++;
		for( int i=0; i<nBin; i++ )
			freqList[i] = freqList[i]/(float)nPoint;
	}

	float
	computeEMD1D( float *dist1, float *dist2, int nBin )
	{
		float dk = 0;
		float emd = 0;
		for(int iB=0; iB<nBin; iB++ )
		{
			//dk += dist1[iB]-dist2[iB];
			dk = dist1[iB]-dist2[iB];
			emd += abs(dk);
		}

		return emd;
	}// End function

	/*
	 *
	 */
	void 
	differentiate( float *data, float *diff, int nP )
	{
		// Using forward difference
		for( int i=0; i<nP-1; i++ )
			diff[i] = data[i+1]-data[i];

	}// End function

	/*
	 *
	 */
	void 
	findMinPoints( float *diff, int nP, int *minP, int *nMinP  )
	{
		int lastSign, curSign;
		*nMinP = 0;

		int minIndex = 0;
		float eps = 0.00001;
		for( int i=0; i<nP; i++ )
		{
			//cout << diff[i] << endl;
			if( i == 0 )
			{
				if( diff[i]>=eps )
					lastSign = 1;
				else if( diff[i]<=eps )
					lastSign = -1;
				else
					lastSign = 0;
			}
			else
			{
				//curSign = ( diff[i]>=0.00001 ? 1 : 0 );
				if( diff[i]>=eps )
					curSign = 1;
				else if( diff[i]<=eps )
					curSign = -1;
				else
					curSign = 0;

				if( curSign == 1 && lastSign == -1 )
				{
					minP[minIndex] = i;
					minIndex++;
					(*nMinP) += 1;
				}
					
				// Update lastsign 
				// for next iteration
				lastSign = curSign;				
			}// End if-else
		}// End for

	}// End function

};
#endif
