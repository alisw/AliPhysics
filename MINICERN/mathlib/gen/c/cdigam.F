*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION CDIGAM(Z)
      COMPLEX CDIGAM,Z
#include "gen/defc64.inc"
     + WDIGAM,W
#include "gen/def64.inc"
     +  D

#include "gen/gcmpfun.inc"
      SROUND(D)=D+(D-SNGL(D))
      W=Z
      W=WDIGAM(W)
      CDIGAM=CMPLX(SROUND(GREAL(W)),SROUND(GIMAG(W)))
      RETURN
      END
#endif
