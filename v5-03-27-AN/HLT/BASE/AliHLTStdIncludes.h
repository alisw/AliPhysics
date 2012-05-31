// @(#) $Id$
// Original file : src/AliHLTStandardIncludes.h,v 1.5

#ifndef ALIHLTSTDINCLUDESH
#define ALIHLTSTDINCLUDESH

#if __GNUC__ >= 3
#include <fstream>
#include <iostream>

#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cerrno>
#include <cstdarg>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT std::cout
#define STDCERR std::cerr
#define STDENDL std::endl
#define STDIF   std::ifstream
#define STDOF   std::ofstream

#else
#include <iostream.h>
#include <fstream.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT cout
#define STDCERR cerr
#define STDENDL endl
#define STDIF   ifstream
#define STDOF   ofstream

#endif //__GNUC__

#endif
