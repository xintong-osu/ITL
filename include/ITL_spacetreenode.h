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
	int level;
	int nDim;
	
	float low[4];
	float high[4];

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
	ITL_spacetreenode( int ndim, int lev, float *l, float *h );

	/**
	 * Destructor
	*/
	~ITL_spacetreenode();

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
	void setEntropy( float ge );

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
	void saveNode( FILE *outFile );
};
#endif
/* ITL_SPACETREENODE_H_ */
