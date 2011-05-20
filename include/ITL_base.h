/**
 *  Base class for the ITL library.
 *  A basic utility class which handles library initialization and related tasks.
 *  Created on: Nov 17, 2010.
 *  @author Abon
 *  @author Teng-Yok
 */
#ifndef ITL_BASE_H
#define ITL_BASE_H

#include "ITL_header.h"

class ITL_base
{
public:

	/**
	 * Library Initialization function.
	 */
	static void ITL_init();
	/**
	 * Function specifies what to do when 'new' operator cannot allocate sufficient memory.
	 */
	static void ITL_new_handler();
};
#endif
/* ITL_BASE_H_ */
