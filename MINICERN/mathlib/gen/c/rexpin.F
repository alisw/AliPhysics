*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:06  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION REXPIN(X)
C     IMPLICIT DOUBLE PRECISION (D)
#include "gen/def64.inc"
     +  D,DEXPIN,DEXPIE

      SROUND(D)=D+(D-SNGL(D))

      ENTRY EXPINT(X)

      REXPIN=SROUND(DEXPIN(DBLE(X)))
      RETURN

      ENTRY REXPIE(X)
      REXPIE=SROUND(DEXPIE(DBLE(X)))
      RETURN
      END
#endif
