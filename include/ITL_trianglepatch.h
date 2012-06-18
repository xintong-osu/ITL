/*
 * ITL_trianglepatch.h
 *
 *  Created on: Apr 18, 2012
 *      Author: abon
 */

#ifndef ITL_TRIANGLEPATCH_H_
#define ITL_TRIANGLEPATCH_H_

#include "ITL_header.h"

class ITL_trianglepatch
{
	int vid[3];
	int nChild;

	//ITL_trianglepatch* child;
	//ITL_trianglepatch* parent;
	int parentID;
	int level;

public:

	ITL_trianglepatch( int a, int b, int c, int l, int p = -1 );
	ITL_trianglepatch( const ITL_trianglepatch& that );
	ITL_trianglepatch& operator= ( const ITL_trianglepatch& that );

	~ITL_trianglepatch();

	int getVertexID( int i );
	int getParentID();

};
#endif
/* ITL_TRIANGLEPATCH_H_ */
