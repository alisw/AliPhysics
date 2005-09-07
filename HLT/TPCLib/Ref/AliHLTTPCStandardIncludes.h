// @(#) $Id$

#ifndef ALIHLTTPCSTANDARDINCLUDESH
#define ALIHLTTPCSTANDARDINCLUDESH

#if __GNUC__>=3
#include <iostream>
#include <fstream>

#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT std::cout
#define STDCERR std::cerr
#define STDENDL std::endl

using namespace std;

#else
#include <iostream.h>
#include <fstream.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

/* Use these only if absolutely necessary 
eg. in inline functions defined in header files */
#define STDCOUT cout
#define STDCERR cerr
#define STDENDL endl

#endif //GCCVERSION

#endif

