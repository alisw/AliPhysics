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
      FUNCTION CLGAMA(Z)
      COMPLEX CLGAMA,Z
      COMPLEX CLOGAM
#include "gen/defc64.inc"
     +     WLGAMA,W
#include "gen/def64.inc"
     +     D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

      SROUND(D)=D+(D-SNGL(D))

      ENTRY CLOGAM(Z)

      W=Z
      W=WLGAMA(W)
      CLGAMA=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      RETURN
      END
#endif
