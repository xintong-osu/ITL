/**
 *  Class for the space-partitioning tree node in ITL library.
 *  A class that implements a node of the space-partitioning data structure. A node represents a partition of the field.
 *  Created on: July 17, 2011.
 *  @author Abon
 *  @author Teng-Yok
 */
#ifndef ITL_SPACETREENODE_H
#define ITL_SPACETREENODE_H

#include "ITL_header.h"

class ITL_spacetreenode
{
	int globalID;

	int level;
	int nDim;
	
	float low[4];
	float high[4];

	int nBin;
	float freqList[720];
	float globalEntropy;
	
	int nChild;
	ITL_spacetreenode *children;
	ITL_spacetreenode *parent;
	
public:

	/**
	 * Default constructor
	 */
	ITL_spacetreenode();

	/**
	 * Constructor
	 */
	ITL_spacetreenode( int nDim );

	/**
	 * Constructor
	 */
	ITL_spacetreenode( int ndim, int lev, float *l, float *h, int nbin = 360 );

	/**
	 * Destructor
	*/
	~ITL_spacetreenode();

	/**
	 * 
	 */
	void setGlobalID( int id );

	/**
	 * 
	 */
	void setLevel( int l );

	/**
	 *
	 */
	void setLimit( float *l, float *h );

	/**
	 * 
	 */
	void setNumDim( int nd );

	/**
	 * 
 	 */
	void setParent( ITL_spacetreenode *p );

	/**
	 *
 	 */
	void setChildren( int nc, ITL_spacetreenode *carray );

	/**
	 *
 	 */
	void setFrequencyList( int nbin, float *freqlist );

	/**
	 *
	 */
	void setEntropy( float ge );

	/**
	 *
	 */
	int getGlobalID();

	/**
	 *
	 */
	int getNumBin();

	/**
	 *
	 */
	int getNumDim();

	/**
	 * 
	 */
	int getLevel();

	/**
	 *
	 */
	int getNumChildren();

	/**
	 * 
 	 */
	ITL_spacetreenode* getParent();

	/**
	 * 
	 */
	ITL_spacetreenode* getChild( int index );

	/**
	 * 
	 */
	float* getLowLimit();

	/**
	 * 
	 */
	float* getHighLimit();

	/**
	 *
 	 */
	float* getFrequencyList();

	/**
	 *
	 */
	float getEntropy();

	/**
	 * 
	 */
	float getVolume();

	/**
	 *
	 */
	void printNode();

	/**
	 *
	 */
	void saveNode( FILE *outFile, int nChildTree );
};
#endif
/* ITL_SPACETREENODE_H_ */
