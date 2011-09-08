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

template <class T>
class ITL_partition
{
	int nDim;
public:
	/**
	 * Default constructor
	 */	
	ITL_partition( int nd )
	{
		nDim = nd;
	}

	/**
	 * Default destructor
	 */	
	~ITL_partition()
	{
	}

	/**
	 *
	*/
	void partition_Octree( float *parentLow, float *parentHigh, ITL_field_regular<T> *parentField, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
	{
		float *mid = new float[nDim];
	
		/*		
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
		*/

		// Compute global entropy of each child 
		computeEntropyChildren( parentField, nchild, nBin, childLow, childHigh, entropyArray );	

		delete [] mid;	
	}// end function

	/**
	 *
	*/
	void partition_Kdtree_SplitMiddle( float *parentLow, float *parentHigh, ITL_field_regular<T> *parentField, int level, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
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
		computeEntropyChildren( parentField, nchild, nBin, childLow, childHigh, entropyArray );	
		
	
	}// end function

	/**
	 *
	*/
	void partition_Kdtree_SplitMiddle2( float *parentLow, float *parentHigh, ITL_field_regular<int> *parentBinField, int level, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
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
		computeEntropyChildren2( parentBinField, nchild, nBin, childLow, childHigh, entropyArray );	
		
	
	}// end function

	/**
	 *
	 */
	void partition_MaxEntropy( float *parentLow, float *parentHigh, ITL_field_regular<T> *parentField, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
	{	
		// Initialize local variables		
		int nTrial = 9;
		//float percntArray[] = {0.2, 0.4, 0.6, 0.8};
		float percntArray[] = {0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66};

		float mid = 0;		
		float entropyArrayCopy[nchild];
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

		printf( "Parent: (%f %f %f) -- (%f %f %f)\n", parentLow[0], parentLow[1], parentLow[2], parentHigh[0], parentHigh[1], parentHigh[2] );
		for( int i = 0; i<nDim-1; i++ )		
		//for( int i = 0; i<nDim; i++ )
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


				printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n", childLow[0][0], childLow[0][1], childLow[0][2], childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n", childLow[1][0], childLow[1][1], childLow[1][2], childHigh[1][0], childHigh[1][1], childHigh[1][2] );

				// Compute global entropy of each child 
				computeEntropyChildren( parentField, nchild, nBin, childLow, childHigh, entropyArray );	
				cout << entropyArray[0] << " " << entropyArray[1] << endl;
				//cout << entropyArrayCopy[0] << " " << entropyArrayCopy[1] << endl;

				// Remember this partition, if it maximizes entropy of the two children				
				//if( max( entropyArray[0], entropyArray[1] ) > max( entropyArrayCopy[0], entropyArrayCopy[1] ) )
				//float score = max( entropyArray[0], entropyArray[1] );
				//float score = ( entropyArray[0]*percntArray[j]+entropyArray[1]*(1-percntArray[j]) );
				//float score = fabs( entropyArray[0]*percntArray[j] - entropyArray[1]*(1-percntArray[j]) );
				float t1 = ( entropyArray[0] / ( (float)log( (float)nBin ) / log(2.0f ) ) ) * percntArray[j];
 				float t2 = ( entropyArray[1] / ( (float)log( (float)nBin ) / log(2.0f ) ) ) * (1-percntArray[j]);
				float score = (t1 + t2)/2.0f + fabs( t1-t2 );
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
	
	}// end function


	/**
	 *
	 */
	void partition_MaxEntropy2( float *parentLow, float *parentHigh, ITL_field_regular<int> *parentBinField, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
	{	
		// Initialize local variables		
		//int nTrial = 4;
		//float percntArray[] = {0.2, 0.4, 0.6, 0.8};
		int nTrial = 9;		
		float percntArray[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
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

				//printf( "Trial child 1: (%f %f %f) -- (%f %f %f)\n", childLow[0][0], childLow[0][1], childLow[0][2], childHigh[0][0], childHigh[0][1], childHigh[0][2] );
 				//printf( "Trial child 2: (%f %f %f) -- (%f %f %f)\n", childLow[1][0], childLow[1][1], childLow[1][2], childHigh[1][0], childHigh[1][1], childHigh[1][2] );

				// Compute global entropy of each child 
				computeEntropyChildren2( parentBinField, nchild, nBin, childLow, childHigh, entropyArray );	
				//cout << entropyArray[0] << " " << entropyArray[1] << endl;
				//cout << entropyArrayCopy[0] << " " << entropyArrayCopy[1] << endl;

				// Remember this partition, if it maximizes entropy of the two children				
				//if( max( entropyArray[0], entropyArray[1] ) > max( entropyArrayCopy[0], entropyArrayCopy[1] ) )
				if( ( entropyArray[0]+entropyArray[1] ) > ( entropyArrayCopy[0]+entropyArrayCopy[1] ) )
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

	void computeEntropyChildren( ITL_field_regular<T> *parentField,  int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
	{
		//int lowPad[3];
		//int highPad[3];
		float lowF[3];
		float highF[3];
		//int sizeNeighborhoodArray[3];
		lowF[0] = lowF[1] = lowF[2] = 0;
		//lowPad[0] = lowPad[1] = lowPad[2] = 0;
		//highPad[0] = highPad[1] = highPad[2] = 0;
		//sizeNeighborhoodArray[0] = sizeNeighborhoodArray[1] = sizeNeighborhoodArray[2] = 0;
		T* childFieldData = NULL;
		ITL_field_regular<T> *childField = NULL;
		ITL_globalentropy<T> *globalEntropyComputer = NULL;
		int freqList[nBin];

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
			childField = new ITL_field_regular<T>( childFieldData, nDim, lowF, highF );//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute global entropy of the sub-field
			globalEntropyComputer = new ITL_globalentropy<T>( childField );
			
			//globalEntropyComputer->computeHistogramBinField( "vector", nBin );			
			globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );

			entropyArray[i] = globalEntropyComputer->getGlobalEntropy();
			globalEntropyComputer->getHistogramFrequencies( nBin, freqList );
			//printf( "%f, %f, %f, %f, %f, %f, ", childLow[i][0], childHigh[i][0], childLow[i][1], childHigh[i][1], childLow[i][2], childHigh[i][2] );   						
			//for( int k=0; k<nBin; k++ )
			//	printf( "%d, ", freqList[k] );
			//printf( "\n" );


			delete globalEntropyComputer;
			delete childField;

		}// end for

		
	}// end function


	void computeEntropyChildren2( ITL_field_regular<int> *parentBinField, int nchild, int nBin, float **childLow, float **childHigh, float *entropyArray )
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

			childField = new ITL_field_regular<T>( nDim, lowF, highF );//, lowPad, highPad, sizeNeighborhoodArray );
			childBinField = new ITL_field_regular<int>( childFieldBinData, nDim, lowF, highF );//, lowPad, highPad, sizeNeighborhoodArray );

			// Compute global entropy of the sub-field
			globalEntropyComputer = new ITL_globalentropy<T>( NULL, childBinField );
			
			globalEntropyComputer->computeHistogramBinField( "vector2", nBin );
			globalEntropyComputer->computeGlobalEntropyOfField( nBin, false, 0 );

			entropyArray[i] = globalEntropyComputer->getGlobalEntropy();

			delete globalEntropyComputer;
			delete childField;
			//delete childBinField;

		}// end for

		
	}// end function

};
#endif
