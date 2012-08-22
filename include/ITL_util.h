/**
 *  Utility class for the ITL library.
 *  A basic utility class which implements basic array operations, timing calculation etc.
 *  Created on: Nov 22, 2010.
 *  @author Abon
 *  @author Teng-Yok
 */

#ifndef ITL_UTIL_H_
#define ITL_UTIL_H_

#include <mpi.h>
#include "ITL_header.h"
#include "ITL_vectormatrix.h"

template <class T>
class ITL_util
{
public:


public:

	static T
	clamp( T val, T low, T high )
	{
		if( val < low)		return low;
		if( val > high )	return high;
		return val;
	}// end function

	static T
	mirror( T val, T low, T high )
	{
		if( val < low)		return low + ( low - val );
		if( val > high )	return high - ( val - high );
		return val;
	}// end function
	
	static T
	Min( T* array, int len )
	{
		T min = array[0];
		for( int i=0; i<len; i++ )
		{
			if( array[i] < min )	min = array[i];
		}
		return min;

	}// end function
	
	static T
	Max( T* array, int len )
	{
		T max = array[0];
		for( int i=0; i<len; i++ )
		{
			if( array[i] > max )	max = array[i];
		}
		return max;

	}// end function

	static T
	sum( T* array, int len )
	{
		T sum = 0;
		//cout << sum << endl;
		for( int i=0; i<len; i++ )
		{
			//cout << array[i] << endl;
			sum += array[i];
		}
		return sum;

	}// end function

	static T
	prod( T* array, int len )
	{
		T product = 1;
		for( int i=0; i<len; i++ )
			product *= array[i];
		return product;

	}// end function

	static void
	fill( T* array, int len, T val )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = val;
		}
	}

	static void
	addArrays( T* one, T* two, T* sum, int len )
	{
		for( int i=0; i<len; i++ )
		{
			sum[i] = one[i]+two[i];
		}
	}

	static void
	subtractArrays( T* one, T* two, T* sub, int len )
	{
		for( int i=0; i<len; i++ )
		{
			sub[i] = one[i]-two[i];
		}
	}

	static void
	addArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]+scalar;
		}
	}

	static void
	subtractArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]-scalar;
		}
	}

	/**
	 * Index conversion: 3D to 1D
	 */
	static int
	index3DTo1D( int x, int y, int z, int* fieldsize )
	{
		return z*fieldsize[1]*fieldsize[0] + y*fieldsize[0] + x;
	}// End function

	/**
	 * Tri-linear interpolation
	 */
	static SCALAR
	triLinterp_scalar( SCALAR array[], float *f )
	{
		SCALAR tmp[2];
		tmp[0] = biLinterp_scalar( array, f );
		tmp[1] = biLinterp_scalar( array+4, f );
		return linterp_scalar( tmp, f[2] );

	}// End function

	/**
	 * Bi-linear interpolation
	 */
	static SCALAR
	biLinterp_scalar( SCALAR array[], float *f )
	{
		SCALAR tmp[2];
		tmp[0] = linterp_scalar( array, f[0] );
		tmp[1] = linterp_scalar( array+1, f[0] );
		return linterp_scalar( tmp, f[1] );

	}// End function

	/**
	 * Linear interpolation
	 */
	static SCALAR
	linterp_scalar( SCALAR array[], float f )
	{
		return ( array[0]*f + array[1]*(1-f) );

	}// End function

	/**
	 * Tri-linear interpolation
	 */
	static VECTOR3
	triLinterp_vector( VECTOR3* array, float *f )
	{
		VECTOR3 tmp[2];
		tmp[0] = biLinterp_vector( array, f );
		tmp[1] = biLinterp_vector( array+4, f );
		return linterp_vector( tmp, f[2] );

	}// End function

	/**
	 * Bi-linear interpolation
	 */
	static VECTOR3
	biLinterp_vector( VECTOR3* array, float *f )
	{
		VECTOR3 tmp[2];

		tmp[0] = linterp_vector( array, f[0] );
		tmp[1] = linterp_vector( array+2, f[0] );

		return linterp_vector( tmp, f[1] );

	}// End function

	/**
	 * Linear interpolation
	 */
	static VECTOR3
	linterp_vector( VECTOR3* array, float f )
	{
		VECTOR3 tmp1( array[0].x()*f, array[0].y()*f, array[0].z()*f );
		VECTOR3 tmp2( array[1].x()*(1-f), array[1].y()*(1-f), array[1].z()*(1-f) );
		VECTOR3 tmp( tmp1[0] + tmp2[0], tmp1[1] + tmp2[1], tmp1[2] + tmp2[2] );

		return tmp;
	}// End function

	/**
	 * Interpolation (nearest neighbor
	 */
	static VECTOR3
	interp_NN_vector( VECTOR3* array, float* f )
	{
		float dist[8];
		dist[0] = sqrt( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] );
		dist[1] = sqrt( (1 - f[0])*(1 - f[0]) + f[1]*f[1] + f[2]*f[2] );
		dist[2] = sqrt( (1 - f[0])*(1 - f[0]) + (1 - f[1])*(1 - f[1]) + f[2]*f[2] );
		dist[3] = sqrt( f[0]*f[0] + (1 - f[1])*(1 - f[1]) + f[2]*f[2] );
		dist[4] = sqrt( f[0]*f[0] + f[1]*f[1] + (1 - f[2])*(1 - f[2]) );
		dist[5] = sqrt( (1 - f[0])*(1 - f[0]) + f[1]*f[1] + (1 - f[2])*(1 - f[2]) );
		dist[6] = sqrt( (1 - f[0])*(1 - f[0]) + (1 - f[1])*(1 - f[1]) + (1 - f[2])*(1 - f[2]) );
		dist[7] = sqrt( f[0]*f[0] + (1 - f[1])*(1 - f[1]) + (1 - f[2])*(1 - f[2]) );

		float mD = dist[0];
		int index = 0;
		for( int i=0; i<8; i++ )
		{
			if( dist[i] < mD )
			{
				mD = dist[i];
				index = i;
			}
		}

		return array[index];

	}// End function

	/**
	 * Start Time computation.
	 */
	static double
	startTimer()
	{
		return MPI_Wtime();
	}// end function

	/**
	 * Complete Time computation.
	 * @param startTime Recored time at start point
	 * @return computed time in seconds from start to end.
	 */
	static double
	endTimer( double startTime )
	{
		double endTime = MPI_Wtime();
		//double resolution = MPI_Wtick();
		return ( endTime - startTime );
	}// end function
	
	static void getArgs( const char *argsFileName, list<string> *argNames, list<string> *argValues )
	{
		char* argName = NULL;// "";
		char* argVal = NULL;//"";
		char nextLine[1000];

		// Open ascii file containing list of arguments
		FILE* argsFile = fopen( argsFileName, "r" );

		// Scan the file to create two lists of strings
		// List of argument names and list of argument values
		while( 1 )
		{
			// Read next line
			fgets( nextLine, 1000, argsFile );

			// Check for EOF
			if( feof( argsFile ) )
				break;

			// Skip to next line if this is a comment or a blank line
			if( *nextLine == '#' || *nextLine == '\n' )
				//printf( "Comment or blank line\n" );
	        		continue;

			// Tokenize string in to argument name and value
			argName = strtok( nextLine, " \n" );
			argVal = strtok( NULL, " \n" );
			string name( argName );
			string val( argVal );
	
			// Push strings in to lists
			argNames->push_back( name );
			argValues->push_back( val );
		}

		// Close file
		fclose( argsFile );

	}// end function

	static char* getArgWithName( const char* name, list<string> *argNames, list<string> *argValues )
	{	
		list<string>::iterator iterName;
		list<string>::iterator iterVal;

		for( iterName = argNames->begin(), iterVal = argValues->begin(); iterName != argNames->end(); ++iterName, ++iterVal )
		{
			if( strcmp( (*iterName).c_str(), name ) == 0 )
				return (char*)( (*iterVal).c_str() );

		}// end for

		return "";

	}// end function

};


#endif
/* ITL_UTIL_H_ */
