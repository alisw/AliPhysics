*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:40  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMHPX)
#include "hpxgs/getarg.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxsys/getarg.F"
#endif
