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
      FUNCTION CWERF(Z)
      COMPLEX CWERF,Z
#include "gen/defc64.inc"
     +     WWERF,W
#include "gen/def64.inc"
     +     D

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

      SROUND(D)=D+(D-SNGL(D))
      W=Z
      W=WWERF(W)
      CWERF=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      END
#endif
