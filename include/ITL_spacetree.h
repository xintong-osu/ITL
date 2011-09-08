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
#define KDTREE_MID 1
#define ENTROPYBASED 2

template <class T>
class ITL_spacetree
{
	int nDim;	
	int nBin;
	
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
	void initTree( int ndim, float *l, float *h, float en = -1 )
	{
		nDim = ndim;	
		nBin = 90;
		root = new ITL_spacetreenode( ndim, 0, l, h );
		if( en != -1 ) root->setEntropy( en );
		curNode = root;
		partitioner = new ITL_partition<T>( ndim );
	}	

	/**
	 * Recursive function that starts from the root and subdivides the spatial domain. 
	 */
	void createTree( ITL_spacetreenode *cur, int partitionType, int level, int maxLevel = 0 )
	{
		
		//cout << "in level " << level << endl;
		// Termination: either by reaching maximum number of levels or by reaching threshold of minimum size
		if( maxLevel != 0 && level == maxLevel ) 
			return;
		if( maxLevel == 0 && cur->getVolume() <= 512.0f )
			return;

		int nChild = -1;
		if( partitionType == OCTREE ) 		  nChild = 8; 		
		else if( partitionType == KDTREE_MID )    nChild = 2;
		else if( partitionType == ENTROPYBASED )  nChild = 2; 		

		// Create partitions based on selected nature
		float **childLow = new float*[nChild];
		float **childHigh = new float*[nChild];
		float *entropyArray = new float[nChild];
		for( int i=0; i<nChild; i++ )
		{
			childLow[i] = new float[nDim];
			childHigh[i] = new float[nDim];
		}
		

		if( partitionType == OCTREE )
		{
			partitioner->partition_Octree( cur->getLowLimit(), cur->getHighLimit(), dataField, nChild, nBin, childLow, childHigh, entropyArray );
		}
		else if( partitionType == KDTREE_MID )
		{
			partitioner->partition_Kdtree_SplitMiddle( cur->getLowLimit(), cur->getHighLimit(), dataField, level, nChild, nBin, childLow, childHigh, entropyArray );
			//partitioner->partition_Kdtree_SplitMiddle2( cur->getLowLimit(), cur->getHighLimit(), binField, level, nChild, nBin, childLow, childHigh, entropyArray );
		}
		else if( partitionType == ENTROPYBASED )
		{
			partitioner->partition_MaxEntropy( cur->getLowLimit(), cur->getHighLimit(), dataField, nChild, nBin, childLow, childHigh, entropyArray );
			//partitioner->partition_MaxEntropy2( cur->getLowLimit(), cur->getHighLimit(), binField, nChild, nBin, childLow, childHigh, entropyArray );
		}

		// Allocate memory for child Nodes
		ITL_spacetreenode *children = new ITL_spacetreenode[nChild];
		
		// Set properties for each child 
		for( int i=0; i<nChild; i++ )
		{
			children[i].setLevel( level+1 );
			children[i].setNumDim( cur->getNumDim() );
			children[i].setEntropy( entropyArray[i] );
			children[i].setLimit( childLow[i], childHigh[i] );
			children[i].setParent( cur );
		}
		
		// Add children to the current node
		cur->setChildren( nChild, children );
		
		// Recursive call to children
		for( int i=0; i<nChild; i++ )
			createTree( cur->getChild(i), partitionType, level+1, maxLevel ); 

		delete [] children;
		delete [] childLow;
		delete [] childHigh;
		delete [] entropyArray;
	}	

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

	/**
	 * 
	 */
	void saveTree( ITL_spacetreenode *cur, FILE *outFile )
	{
		// Print current node
		cur->saveNode( outFile );

		// Termination
		if( cur->getNumChildren() == 0 )
			return;

		// Recursive call to children
		for( int i=0; i<cur->getNumChildren(); i++ )
			saveTree( cur->getChild(i), outFile );
	}// end function		
	
};
#endif
/* ITL_SPACETREE_H_ */
