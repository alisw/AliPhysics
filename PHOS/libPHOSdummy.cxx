#include <stdio.h>
/*

    dummy entry points for phos

    */

#ifdef WIN32
#define shlist SHLIST
#define shinit SHINIT
#define shevnt SHEVNT
#define reconsfirst RECONSFIRST
#define type_of_call  _stdcall
#else
#define shlist shlist_
#define shinit shinit_
#define shevnt shevnt_
#define reconsfirst reconsfirst_
#define type_of_call
#endif

#define DUMMY(name) \
extern "C" type_of_call void name() {\
     printf("Dummy version of \"" #name "\" reached \n"); \
     }

DUMMY(shlist)
DUMMY(shinit)
DUMMY(shevnt)
DUMMY(reconsfirst)


