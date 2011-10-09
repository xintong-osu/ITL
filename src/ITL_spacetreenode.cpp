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
ITL_spacetreenode::ITL_spacetreenode( int ndim, int lev,
									  float *l, float *h,
									  int nbin )
{
	nDim = ndim;
	nBin = nbin;
	assert( nBin <= 720 );
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
}// End destructor

/**
 * 
 */
void
ITL_spacetreenode::setGlobalID( int id )
{
	globalID = id;
}// End function

/**
 * 
 */
void
ITL_spacetreenode::setLimit( float *l, float *h )
{
	memcpy( low, l, nDim * sizeof(float) );
	memcpy( high, h, nDim * sizeof(float) ); 

}// End function

/**
 * 
 */
void
ITL_spacetreenode::setLevel( int l )
{
	level = l;
}// End function

/**
 *
 */
void
ITL_spacetreenode::setNumDim( int nd )
{
	nDim = nd;
}// End function

/**
 * 
 */
void
ITL_spacetreenode::setParent( ITL_spacetreenode *p )
{
	parent = p;
}

/**
 * 
 */
void
ITL_spacetreenode::setChildren( int nc, ITL_spacetreenode *carray )
{
	nChild = nc;
	
	if( children != NULL ) delete [] children;
	children = new ITL_spacetreenode[nChild];

	//memcpy( children, carray, sizeof( ITL_spacetreenode ) * nChild );	 
	// Set properties for each child 
	for( int i=0; i<nChild; i++ )
	{
		children[i].setLevel( carray[i].getLevel() );
		children[i].setNumDim( carray[i].getNumDim() );
		children[i].setFrequencyList( nBin, carray[i].getFrequencyList() );
		children[i].setEntropy( carray[i].getEntropy() );
		children[i].setLimit( carray[i].getLowLimit(), carray[i].getHighLimit() );
		children[i].setParent( this );
	}

}// End function

/**
 *
 */
void
ITL_spacetreenode::setFrequencyList( int nbin, float *freqlist )
{
	assert( nBin <= 720 );
	nBin = nbin;
	memcpy( this->freqList, freqlist, sizeof(float)*nBin );
	
}// End function

/**
 *
 */
void
ITL_spacetreenode::setEntropy( float ge )
{
	globalEntropy = ge;
}// End function

/**
 * 
 */
int
ITL_spacetreenode::getGlobalID()
{
	return globalID;
}// End function

/**
 * 
 */
int
ITL_spacetreenode::getLevel()
{
	return level;
}

/**
 * 
 */
int
ITL_spacetreenode::getNumDim()
{
	return nDim;
}

/**
 * 
 */
int
ITL_spacetreenode::getNumBin()
{
	return nBin;
}

/**
 * 
 */
int
ITL_spacetreenode::getNumChildren()
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
float*
ITL_spacetreenode::getFrequencyList()
{
	return this->freqList;
}// end function

/**
 *
 */
float
ITL_spacetreenode::getEntropy()
{
	return globalEntropy;
}// end function

/**
 * 
 */
float
ITL_spacetreenode::getVolume()
{
	float volume = 1;
	for( int i=0; i<nDim; i++ )
		volume *= ( high[i] - low[i] );

	return volume;
}// end function

/**
 *
 */
void
ITL_spacetreenode::printNode()
{
	printf( "%d %d (%f %f) (%f %f) (%f %f), %f\n",
			level, globalID,
			low[0], high[0], low[1],
			high[1], low[2], high[2],
			globalEntropy );
}// end function

/**
 *
 */
void
ITL_spacetreenode::saveNode( FILE *outFile, int nChildTree )
{
	// Global ID
	fprintf( outFile, "%d, ", globalID );

	// Level
	fprintf( outFile, "%d, ", level );

	// Pointer to parent
	if( this->parent == NULL )
		fprintf( outFile, "-1, " );
	else
		fprintf( outFile, "%d, ", this->parent->getGlobalID() );

	// Pointers to children
	for( int i=0; i<nChildTree; i++ )
	{
		if( this->children == NULL )
			fprintf( outFile, "-1, " );
		else
			fprintf( outFile, "%d, ", this->children[i].getGlobalID() );
	}
	
	// Geometric extent
	fprintf( outFile, "%f, %f, %f, %f, %f, %f, ",
			 low[0], high[0], low[1],
			 high[1], low[2], high[2] );

	// Entropy
	fprintf( outFile, "%f, ", getEntropy() );

	// Histogram
	for( int i=0; i<nBin-1; i++ )
		fprintf( outFile, "%f, ", freqList[i] );
	fprintf( outFile, "%f\n", freqList[nBin-1] );

}// end function






