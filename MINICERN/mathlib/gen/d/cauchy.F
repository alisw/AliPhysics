*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CAUCHY(F,A,B,S,EPS)
      CAUCHY=RCAUCH(F,A,B,S,EPS)
      RETURN
      END
#endif
