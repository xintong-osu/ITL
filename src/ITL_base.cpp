/**
 * @file ITL_base.cpp
 * Source file for ITL_base.
 * Created on: Nov 18, 2010
 * @author Abon
 * @author Teng-Yok
 */

#include "ITL_base.h"

/**
 * Library Initialization function.
 */
void ITL_base::ITL_init()
{
	// set new handler
	set_new_handler( ITL_base::ITL_new_handler );

}// end function

/**
 * Function specifies what to do when 'new' operator cannot allocate sufficient memory.
 */
void ITL_base::ITL_new_handler()
{
  printf( "Error: Failed to allocate memory\n" );
  exit (1);

}// end function
