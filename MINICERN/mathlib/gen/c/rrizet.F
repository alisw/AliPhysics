*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RRIZET(X)
      DOUBLE PRECISION D,DRIZET

      SROUND(D)=D+(D-SNGL(D))
      RRIZET=SROUND(DRIZET(DBLE(X)))
      RETURN
      END
#endif
