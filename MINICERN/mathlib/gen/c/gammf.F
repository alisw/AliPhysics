*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION GAMMF(X)
      DOUBLE PRECISION DGAMMF,D

      SROUND(D)=D+(D-SNGL(D))
      GAMMF=SROUND(DGAMMF(DBLE(X)))
      RETURN
      END
#endif
