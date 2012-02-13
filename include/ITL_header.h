/*
 * A Collection of all includes and definitions.
 * Adopted from OSU Flow Vis Library
 * Created: Han-Wei Shen, Liya Li
 * The Ohio State University
 * Date: 06/2005
 *
 */

#ifndef _ITL_HEADER_H_
#define _ITL_HEADER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#if	0	// DEL-BY-LEETEN 02/13/2012-BEGIN
	#if defined( _WIN32 ) || defined( _WIN64 )
		#define _USE_MATH_DEFINES
	#endif
#endif	// DEL-BY-LEETEN 02/13/2012-END
#include <cmath>
#include <cassert>
#include <float.h>
#include <ctime>
#include <vector>
#include <string>
#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <iterator>

using namespace std;

const double	pi = 3.14159265358979323846;			// PI
const double	eps = 1.0E-6;
//const double	DEG_TO_RAD = 0.0174532925199432957692;		// PI / 180
//const double	RAD_TO_DEG = 57.2957795130823208768;		// 180 / PI
//const double	PIBY2 = 1.57079632679489661923;			// PI / 2
//const int	OCT = 8;

#define CLOCKS_PER_MS CLOCKS_PER_SEC/1000.0f
//#define DEBUG_MODE

//#define DUPLICATE
#define MIRROR

#define VERY_HIGH_VALUE 1.0E6;

#endif
