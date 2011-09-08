#include "ITL_spacetreenode.h"

/**
 * Default constructor
 */
ITL_spacetreenode::ITL_spacetreenode()
{
	nDim = 3;
	nChild = 0;
	children = NULL;
	parent = NULL;

}// end constructor

/**
 * Constructor 2
 */
ITL_spacetreenode::ITL_spacetreenode( int ndim )
{
	nDim = ndim;
	nChild = 0;
	children = NULL;
	parent = NULL;

}// end constructor

/**
 * Constructor 3
 */
ITL_spacetreenode::ITL_spacetreenode( int ndim, int lev, float *l, float *h )
{
	nDim = ndim;
	nChild = 0;
	level = lev;
	memcpy( low, l, nDim * sizeof(float) );
	memcpy( high, h, nDim * sizeof(float) ); 

	children = NULL;
	parent = NULL;

}// end constructor

/**
 * Destructor
 */
ITL_spacetreenode::~ITL_spacetreenode()
{
	//if( parent != NULL ) delete parent;	
	for( int i=0; i<nChild; i++ )
		delete ( children+i );		
}// end destructor

/**
 * 
 */
void ITL_spacetreenode::setLimit( float *l, float *h )
{
	memcpy( low, l, nDim * sizeof(float) );
	memcpy( high, h, nDim * sizeof(float) ); 

}// end function

/**
 * 
 */
void ITL_spacetreenode::setLevel( int l )
{
	level = l;
}

/**
 *
 */
void ITL_spacetreenode::setNumDim( int nd )
{
	nDim = nd;
}// end function

/**
 * 
 */
void ITL_spacetreenode::setParent( ITL_spacetreenode *p )
{
	parent = p;
}

/**
 * 
 */
void ITL_spacetreenode::setChildren( int nc, ITL_spacetreenode *carray )
{
	nChild = nc;
	
	if( children != NULL ) delete [] children;
	children = new ITL_spacetreenode[nChild];

	memcpy( children, carray, sizeof( ITL_spacetreenode ) * nChild );	 

}// end function

/**
 *
 */
void ITL_spacetreenode::setEntropy( float ge )
{
	globalEntropy = ge;
}// end function

/**
 * 
 */
int ITL_spacetreenode::getLevel()
{
	return level;
}

/**
 * 
 */
int ITL_spacetreenode::getNumDim()
{
	return nDim;
}

/**
 * 
 */
int ITL_spacetreenode::getNumChildren()
{
	return nChild;
}

/**
 *
 */
float* ITL_spacetreenode::getLowLimit()
{
	return &low[0];
}// end function

/**
 * 
 */
float* ITL_spacetreenode::getHighLimit()
{
	return &high[0];
}// end function

/**
 * 
 */
ITL_spacetreenode* ITL_spacetreenode::getParent()
{
	return parent;
}// end function

/**
 *
 */
ITL_spacetreenode* ITL_spacetreenode::getChild( int index )
{
	return (children+index);
}// end function

/**
 *
 */
float ITL_spacetreenode::getEntropy()
{
	return globalEntropy;
}// end function

/**
 * 
 */
float ITL_spacetreenode::getVolume()
{
	float volume = 1;
	for( int i=0; i<nDim; i++ )
		volume *= ( high[i] - low[i] );

	return volume;
}// end function

/**
 *
 */
void ITL_spacetreenode::printNode()
{
	printf( "%d: (%f %f) (%f %f) (%f %f), %f\n", level, low[0], high[0], low[1], high[1], low[2], high[2], globalEntropy );
}// end function

/**
 *
 */
void ITL_spacetreenode::saveNode( FILE *outFile )
{
	fprintf( outFile, "%f, %f, %f, %f, %f, %f, %d, %f\n", low[0], high[0], low[1], high[1], low[2], high[2], level, globalEntropy );
}// end function






