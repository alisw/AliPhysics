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
      FUNCTION CGPLG(N,M,X)
C
      COMPLEX    CGPLG
#include "gen/defc64.inc"
     +   WGPLG,W
#include "gen/def64.inc"
     +   D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
#include "gen/gcmpfun.inc"
      SROUND(D)=D+(D-SNGL(D))
      W=WGPLG(N,M,DBLE(X))
      CGPLG =CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      RETURN
      END
#endif
