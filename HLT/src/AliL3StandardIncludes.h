#ifndef ALIL3STANDARDINCLUDESH
#define ALIL3STANDARDINCLUDESH

#if GCCVERSION == 3

#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>

#include <iostream>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT std::cout
#define STDCERR std::cerr
#define STDENDL std::endl

#else

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include <iostream.h>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT cout
#define STDCERR cerr
#define STDENDL endl

#endif //GCCVERSION

#endif

