*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (defined(CERNLIB_DOUBLE))&&(!defined(CERNLIB_IBM))
      FUNCTION GAMMA(X)
      DOUBLE PRECISION DGAMMA,D

      SROUND(D)=D+(D-SNGL(D))
      GAMMA=SROUND(DGAMMA(DBLE(X)))
      RETURN
      END
#endif
