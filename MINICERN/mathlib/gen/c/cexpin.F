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
      FUNCTION CEXPIN(Z)
      COMPLEX CEXPIN,Z
#include "gen/defc64.inc"
     +     WEXPIN,W
#include "gen/def64.inc"
     +     D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

      SROUND(D)=D+(D-SNGL(D))
      W=Z
      W=WEXPIN(W)
      CEXPIN=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      RETURN
      END
#endif
