/**
 *  Class for the space-partitioning tree in ITL library.
 *  A class that implements a basic space-partitioning data structure.
 *  Created on: July 17, 2011.
 *  @author Abon
 *  @author Teng-Yok
 */
#ifndef ITL_SPACETREE_H
#define ITL_SPACETREE_H

#include "ITL_header.h"
#include "ITL_spacetreenode.h"
#include "ITL_field_regular.h"

#define OCTREE 0
#define QUADTREE 1
#define KDTREE_MID 2
#define KDTREE_MID_2D 3
#define ENTROPYBASED 4
#define EMDBASED 5

template <class T>
class ITL_spacetree
{
	int nDim;
	int nBin;
	int numChild;

	ITL_partition<T> *partitioner;
	ITL_spacetreenode *root;
	ITL_spacetreenode *curNode;

	ITL_field_regular<T> *dataField;
	ITL_field_regular<int> *binField;

public:

	/**
	 * Constructor
	 */
	ITL_spacetree()
	{
		partitioner = NULL;
		root = NULL;
		curNode = NULL;
		dataField = NULL;
		binField = NULL;
	}

	/**
	 * Destructor
	 */
	~ITL_spacetree()
	{
	}

	/**
	 * Initialization of tree.
	 * Create a single node tree that contains the entire spatial domain
	 */
	void initTree( int ndim, int mappingType,
				   int nbin,
				   float *l, float *h,
				   float *rootfreqlist,
				   float en = -1 )
	{
		nDim = ndim;
		nBin = nbin;
		root = new ITL_spacetreenode( ndim, 0, l, h );
		if( en != -1 ) root->setEntropy( en );
		root->setFrequencyList( nBin, rootfreqlist );
		curNode = root;
		partitioner = new ITL_partition<T>( ndim, mappingType );
	}

	/**
	 * Recursive function that starts from the root and
	 * subdivides the spatial domain.
	 */
	void createTree( ITL_spacetreenode *cur, int partitionType,
					 int level, int maxLevel = 0,
					 FILE *dumpFile = NULL )
	{

		//cout << "in level " << level << endl;
		// Termination: either by reaching maximum number of
		// levels or by reaching threshold of minimum size
		if( maxLevel != 0 && level == maxLevel )
			return;
		//if( maxLevel == 0 && cur->getVolume() <= 512.0f )
		//	return;

		if( partitionType == OCTREE ) 				numChild = 8;
		if( partitionType == QUADTREE ) 			numChild = 4;
		else if( partitionType == KDTREE_MID )		numChild = 2;
		else if( partitionType == KDTREE_MID_2D )	numChild = 2;
		else if( partitionType == ENTROPYBASED )	numChild = 2;
		else if( partitionType == EMDBASED )		numChild = 2;

		// Create partitions based on selected nature
		float **childLow = new float*[numChild];
		float **childHigh = new float*[numChild];
		float **childFreqList = new float*[numChild];
		float *entropyArray = new float[numChild];
		for( int i=0; i<numChild; i++ )
		{
			childLow[i] = new float[nDim];
			childHigh[i] = new float[nDim];
			childFreqList[i] = new float[nBin];
		}

		int stopPartition = 0;
		if( partitionType == OCTREE )
		{
			partitioner->partition_Octree( cur->getLowLimit(), cur->getHighLimit(),
										   dataField, numChild, nBin,
										   childLow, childHigh,
										   childFreqList, entropyArray );
		}
		else if( partitionType == QUADTREE )
		{
			partitioner->partition_Quadtree( cur->getLowLimit(), cur->getHighLimit(),
											 dataField, numChild, nBin,
											 childLow, childHigh,
											 childFreqList, entropyArray );
		}
		else if( partitionType == KDTREE_MID )
		{
			partitioner->partition_Kdtree_SplitMiddle( cur->getLowLimit(), cur->getHighLimit(),
													   dataField, level, 3,
													   numChild, nBin,
													   childLow, childHigh,
													   childFreqList, entropyArray );
			//partitioner->partition_Kdtree_SplitMiddle2( cur->getLowLimit(), cur->getHighLimit(),
			//											binField, level, numChild, nBin,
			//											childLow, childHigh,
			//											childFreqList, entropyArray );
		}
		else if( partitionType == KDTREE_MID_2D )
		{
			//partitioner->partition_Kdtree_SplitMiddle( cur->getLowLimit(), cur->getHighLimit(),
			//										   dataField, level, 2,
			//										   numChild, nBin,
			//										   childLow, childHigh,
			//										   childFreqList, entropyArray );
			partitioner->partition_Kdtree_SplitMiddle2( cur->getLowLimit(), cur->getHighLimit(),
														binField, level, 2,
														numChild, nBin,
														childLow, childHigh,
														childFreqList, entropyArray );
		}
		else if( partitionType == ENTROPYBASED )
		{
			partitioner->partition_MaxEntropy( cur->getLowLimit(), cur->getHighLimit(),
											   dataField, numChild, nBin,
											   childLow, childHigh,
											   childFreqList, entropyArray,
											   dumpFile );
			//partitioner->partition_MaxEntropy2( cur->getLowLimit(), cur->getHighLimit(),
			//								      binField, nChild, nBin,
			//									  childLow, childHigh,
			//									  childFreqList, entropyArray );
		}
		else if( partitionType == EMDBASED )
		{
			//stopPartition = partitioner->partition_EMD(
			//							cur->getLowLimit(),
			//							cur->getHighLimit(),
			//						    dataField, numChild, nBin,
			//							childLow, childHigh,
			//							childFreqList, entropyArray,
			//							dumpFile );
			stopPartition = partitioner->partition_EMD2(
										cur->getLowLimit(),
										cur->getHighLimit(),
									    binField, numChild, nBin,
										childLow, childHigh,
										childFreqList, entropyArray,
										dumpFile );
		}// End if-else

		ITL_spacetreenode *children = NULL;
		if( stopPartition == 0 )
		{

			// Allocate memory for child Nodes
			children = new ITL_spacetreenode[numChild];

			// Set properties for each child
			for( int i=0; i<numChild; i++ )
			{
				children[i].setLevel( level+1 );
				children[i].setNumDim( cur->getNumDim() );
				children[i].setFrequencyList( nBin, childFreqList[i] );
				children[i].setEntropy( entropyArray[i] );
				children[i].setLimit( childLow[i], childHigh[i] );
				children[i].setParent( cur );
			}

			// Add children to the current node
			cur->setChildren( numChild, children );

			// Recursive call to children
			for( int i=0; i<numChild; i++ )
				createTree( cur->getChild(i), partitionType, level+1, maxLevel, dumpFile );
		}
		else
			cur->setChildren( 0, NULL );


		if( children != NULL )
			delete [] children;
		delete [] childLow;
		delete [] childHigh;
		delete [] childFreqList;
		delete [] entropyArray;

	}// End children

	/**
	 *
	 */
	void addNode()
	{
	}// end function

	/**
	 *
	 */
	void deleteNode()
	{
	}// end function

	/**
	 *
	 */
	void setDataField( ITL_field_regular<T> *datafield  )
	{
		dataField = datafield;
	}// end function

	/**
	 *
	 */
	void setBinField( ITL_field_regular<int> *binfield  )
	{
		binField = binfield;
	}// end function

	/**
	 *
	 */
	void setEntropyParameters( int nb  )
	{
		nBin = nb;
	}// end function

	/**
	 *
	 */
	ITL_spacetreenode* getRoot()
	{
		return root;
	}// end function

	/**
	 *
	 */
	ITL_spacetreenode* getCurNode()
	{
		return curNode;
	}// end function

	/**
	 *
	 */
	int setNodeIDTree( ITL_spacetreenode *cur, int latestID )
	{
		// set ID of the current node
		cur->setGlobalID( latestID );

		// Termination
		if( cur->getNumChildren() == 0 )
			return (latestID+1);

		// Recursive call to children
		// Each one increments the global ID
		// and passes back the value after
		// incrementing by 1.
		int nextAvailableID = latestID+1;
		for( int i=0; i<cur->getNumChildren(); i++ )
		{
			nextAvailableID = setNodeIDTree( cur->getChild(i), nextAvailableID );
		}

		return nextAvailableID;
	}// end function

	/**
	 *
	 */
	void saveLeafNodes( ITL_spacetreenode *cur, FILE *outFile )
	{
		// Print current node
		if( cur->getNumChildren() == 0 )
			cur->saveNode( outFile, numChild );

		// Termination
		if( cur->getNumChildren() == 0 )
			return;

		// Recursive call to children
		for( int i=0; i<cur->getNumChildren(); i++ )
			saveLeafNodes( cur->getChild(i), outFile );
	}// end function

	/**
	 *
	 */
	void saveTree( ITL_spacetreenode *cur, FILE *outFile )
	{
		// Print current node
		cur->saveNode( outFile, numChild );

		// Termination
		if( cur->getNumChildren() == 0 )
			return;

		// Recursive call to children
		for( int i=0; i<cur->getNumChildren(); i++ )
			saveTree( cur->getChild(i), outFile );
	}// end function

	/**
	 *
	 */
	void printTree( ITL_spacetreenode *cur )
	{
		// Print current node
		cur->printNode();

		// Termination
		if( cur->getNumChildren() == 0 )
			return;

		// Recursive call to children
		for( int i=0; i<cur->getNumChildren(); i++ )
			printTree( cur->getChild(i) );
	}// end function


};
#endif
/* ITL_SPACETREE_H_ */
