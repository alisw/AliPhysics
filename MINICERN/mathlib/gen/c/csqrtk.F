*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      COMPLEX FUNCTION CSQRTK(Z,N,M)
      COMPLEX Z,SZ
      NF=N+5
      SZ=SQRT(Z)
      S=ABS(AIMAG(SZ))
      GOTO (21,22),M
   21 GOTO (9,2,2,9,11,11,2,2,11),NF
   22 GOTO (2,11,2,2,2,2,2,9,2),NF
    9 SZ=CMPLX(0.0,-S)
      GOTO 2
   11 SZ=CMPLX(0.0,S)
    2 CSQRTK=SZ
      RETURN
      END
