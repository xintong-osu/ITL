/**
 *  Statistical utility class for the ITL library.
 *  A basic utility class which implements basic statistical operations such as mean and standard deviation.
 *  Created on: May 03, 2011.
 *  @author Abon
 *  @author Teng-Yok
 */

#ifndef ITL_STATUTIL_H_
#define ITL_STATUTIL_H_

#include "ITL_header.h"
#include "ITL_util.h"

template <class T>
class ITL_statutil
{
public:

	static T
	Mean( T* array, int len )
	{
		T sum = ITL_util<T>::sum( array, len );
		//cout << sum << endl;
		return ( sum / (double)len );
		
	}// end function

	static T
	Variance( T* array, int len, T mean )
	{
		T mse = 0;
		for( int i=0; i<len; i++ )
			mse += ( array[i] - mean ) * ( array[i] - mean ) ;
		
		return ( mse / (float)len );

	}// end function

	static T
	STD( T* array, int len, T mean )
	{
		return sqrt( ITL_statutil<T>::Variance( array, len, mean ) );

	}// end function	

};


#endif 
/* ITL_STATUTIL_H_ */
