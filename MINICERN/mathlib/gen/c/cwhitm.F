*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION CWHITM(Z,KA,MU)
      COMPLEX CWHITM,Z,KA,MU
#include "gen/defc64.inc"
     + WWHITM,WZ,WKA,WMU,W
#include "gen/def64.inc"
     + D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
#include "gen/gcmpfun.inc"
      SROUND(D)=D+(D-SNGL(D))
      WZ=Z
      WKA=KA
      WMU=MU
      W=WWHITM(WZ,WKA,WMU)
      CWHITM=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      RETURN
      END
#endif
