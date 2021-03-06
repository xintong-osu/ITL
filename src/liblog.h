#ifndef __LIB_LOG_H__
#define __LIB_LOG_H__

// ADD-BY-LEETEN 12/30/2009-BEGIN
#include <iostream>
using namespace std;
// ADD-BY-LEETEN 12/30/2009-END


#include <stdio.h>
#include <stdarg.h>

#ifdef WIN32
	#include <string.h>
	// MOD-BY-LEETEN 02/16/2010-FROM:
		// #define PRINT_HEADER	printf("[%s() @ %s(%d)]: ", __FUNCTION__, strrchr(__FILE__, '\\') + 1, __LINE__)	
	// TO:
	#define PRINT_HEADER	cout<<"["<<__FUNCTION__<<"() @ "<<strrchr(__FILE__, '\\') + 1<<"("<<__LINE__<<")]: ";
	// MOD-BY-LEETEN 02/16/2010-END
#else
	// MOD-BY-LEETEN 02/16/2010-FROM:
		// #define PRINT_HEADER	printf("[%s()@%s(%d)]: ", __FUNCTION__, __FILE__, __LINE__)	
	// TO:
	#define PRINT_HEADER	cout<<"["<<__FUNCTION__<<"() @ "<<__FILE__<<"("<<__LINE__<<")]: ";
	// MOD-BY-LEETEN 02/16/2010-END
#endif

#define	LOG_VAR(var)					\
	{									\
		PRINT_HEADER;					\
		cout<<#var<<" = "<<(var)<<endl;	\
	}									\

#if	0	// MOD-BY-LEETEN 02/16/2010-FROM:
	#define LOG(printf_stmt)	\
		{						\
			PRINT_HEADER;		\
			(printf_stmt);		\
			printf("\n");		\
		}						\

#else	// MOD-BY-LEETEN 02/16/2010-TO:

#define LOG(printf_stmt)	\
	{						\
		PRINT_HEADER;		\
		(printf_stmt);		\
		cout<<endl;			\
	}						

#endif	// MOD-BY-LEETEN 02/16/2010-END

// ADD-BY-LEETEN 03/22/2010-BEGIN

// ADD-BY-LEETEN 06/28/2010-BEGIN
#ifdef WIN32	
// ADD-BY-LEETEN 06/28/2010-END
#define PRINT_ERROR_HEADER	cerr<<"["<<__FUNCTION__<<"() @ "<<strrchr(__FILE__, '\\') + 1<<"("<<__LINE__<<")]: ";
// ADD-BY-LEETEN 06/28/2010-BEGIN
#else
	#define PRINT_ERROR_HEADER	cerr<<"["<<__FUNCTION__<<"() @ "<<__FILE__<<"("<<__LINE__<<")]: ";
#endif
// ADD-BY-LEETEN 06/28/2010-END

#define LOG_ERROR(printf_stmt)	\
	{						\
		PRINT_ERROR_HEADER;		\
		(printf_stmt);		\
		cerr<<endl;			\
	}						

#define LOG_VAR_TO_ERROR(var)	\
	{						\
		PRINT_ERROR_HEADER;		\
		cerr<<#var<<" = "<<(var)<<endl;	\
	}						

// ADD-BY-LEETEN 03/22/2010-END

// ADD-BY-LEETEN 06/28/2010-BEGIN
#define ASSERT_OR_LOG(condition, printf_stmt)		\
	{						\
		if( ! (condition) )			\
		{					\
			char *ASSERT = #condition;	\
			LOG_VAR_TO_ERROR(ASSERT);	\
			LOG_ERROR(printf_stmt);		\
			exit(EXIT_FAILURE);		\
		}					\
	}						\

// ADD-BY-LEETEN 06/28/2010-END

#endif	// __LIB_LOG_H__

/*

$Log: liblog.h,v $
Revision 1.5  2011-07-06 03:22:46  leeten

[07/05/2011]
1. [MOD] solve the bug that crashes LOG_ERROR and LOG_VAR_TO_ERROR.
2. [ADD] Define a new macro ASSERT_OR_LOG.

Revision 1.4  2010-03-25 15:43:03  leeten

[03/23/2010]
1. [ADD] Define new macro LOG_ERROR and LOG_VAR_TO_ERROR to output the log to standard error.

Revision 1.3  2010/02/16 20:03:51  leeten

[02/16/2010]
1. [MOD] Use cout instead printf.

Revision 1.2  2009/12/31 01:45:35  leeten

[12/30/2009]
1. [ADD] Include the header iostream and use the namespace std.

Revision 1.1  2009/06/03 21:56:54  leeten

[2009/06/03]
1. [1ST] This project define a library to log out statement and variables.


*/
