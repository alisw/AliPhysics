*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RDILOG(X)
      DOUBLE PRECISION D,DDILOG

      SROUND(D)=D+(D-SNGL(D))

      RDILOG=SROUND(DDILOG(DBLE(X)))
      RETURN

      ENTRY DILOG(X)

      DILOG=SROUND(DDILOG(DBLE(X)))
      RETURN

      END
#endif
