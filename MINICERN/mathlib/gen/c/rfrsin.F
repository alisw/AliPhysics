*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:02  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RFRSIN(X)
#include "gen/def64.inc"
     +    DFRSIN,DFRCOS,D

      SROUND(D)=D+(D-SNGL(D))

      RFRSIN=SROUND(DFRSIN(DBLE(X)))
      RETURN

      ENTRY FRSIN(X)
      FRSIN=SROUND(DFRSIN(DBLE(X)))
      RETURN

      ENTRY RFRCOS(X)
      RFRCOS=SROUND(DFRCOS(DBLE(X)))
      RETURN

      ENTRY FRCOS(X)
      FRCOS=SROUND(DFRCOS(DBLE(X)))
      RETURN

      END
#endif
