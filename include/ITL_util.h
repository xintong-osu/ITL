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

	static T clamp( T val, T low, T high )
	{
		if( val < low)		return low;
		if( val > high )	return high;
		return val;
	}// end function

	static T mirror( T val, T low, T high )
	{
		if( val < low)		return low + ( low - val );
		if( val > high )	return high - ( val - high );
		return val;
	}// end function
	
	static T Min( T* array, int len )
	{
		T min = array[0];
		for( int i=0; i<len; i++ )
		{
			if( array[i] < min )	min = array[i];
		}
		return min;

	}// end function
	
	static T Max( T* array, int len )
	{
		T max = array[0];
		for( int i=0; i<len; i++ )
		{
			if( array[i] > max )	max = array[i];
		}
		return max;

	}// end function


	static T sum( T* array, int len )
	{
		T sum = 0;
		for( int i=0; i<len; i++ )
			sum += array[i];
		return sum;

	}// end function

	static T prod( T* array, int len )
	{
		T product = 1;
		for( int i=0; i<len; i++ )
			product *= array[i];
		return product;

	}// end function

	static void fill( T* array, int len, T val )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = val;
		}
	}

	static T* addArrays( T* one, T* two, int len )
	{
		T* sum = new T[len];
		for( int i=0; i<len; i++ )
		{
			sum[i] = one[i]+two[i];
		}
		return sum;
	}

	static T* subtractArrays( T* one, T* two, int len )
	{
		T* sub = new T[len];
		for( int i=0; i<len; i++ )
		{
			sub[i] = one[i]-two[i];
		}
		return sub;
	}

	static void addArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]+scalar;
		}
	}

	static void subtractArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]-scalar;
		}
	}

	/**
	 * Start Time computation.
	 */
	static double startTimer()
	{
		return MPI_Wtime();
	}// end function

	/**
	 * Complete Time computation.
	 * @param startTime Recored time at start point
	 * @return computed time in seconds from start to end.
	 */
	static double endTimer( clock_t startTime )
	{
		double endTime = MPI_Wtime();
		//double resolution = MPI_Wtick();
		return ( endTime - startTime );// /  CLOCKS_PER_MS;
	}// end function


};


#endif
/* ITL_UTIL_H_ */
