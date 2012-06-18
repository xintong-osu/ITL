/*
 * ITL_trianglepatch.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: abon
 */

#include "ITL_trianglepatch.h"

ITL_trianglepatch::ITL_trianglepatch( int a, int b, int c, int l, int p )
{
	vid[0] = a;
	vid[1] = b;
	vid[2] = c;
	level = l;
	parentID = p;
	nChild = 0;

	//child = NULL;
	//parent = NULL;

}// end constructor

ITL_trianglepatch::ITL_trianglepatch( const ITL_trianglepatch& that )
{
	for( int i=0; i<3; i++ )
		this->vid[i] = that.vid[i];
	this->level = that.level;
	this->nChild = that.nChild;
	this->parentID = that.parentID;

	//this->child = NULL;
	//this->parent = NULL;

}// end copy constructor

ITL_trianglepatch&
ITL_trianglepatch::operator= ( const ITL_trianglepatch& that )
{
	if ( this != &that ) // protect against invalid self-assignment
	{
		for( int i=0; i<3; i++ )
			vid[i] = that.vid[i];
		level = that.level;
		nChild = that.nChild;
		parentID = that.parentID;

		//child = NULL;
		//parent = NULL;
	}
	// by convention, always return *this
	return *this;
}

ITL_trianglepatch::~ITL_trianglepatch()
{
	//if( child != NULL )
	//	delete [] child;
}

int
ITL_trianglepatch::getVertexID( int i )
{
	return vid[i];
}

int
ITL_trianglepatch::getParentID()
{
	return parentID;
}

