*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      COMPLEX FUNCTION CLOGOK(Z,N,M)
      COMPLEX Z,SZ
      DATA PI/3.1415926535898/
      SZ=LOG(Z)
      S=REAL(SZ)
      NF=N+5
      GOTO (21,22,23) ,M
   21 GOTO (2,2,2,11,10,9,2,2,2),NF
   22 GOTO (2,11,2,2,2,2,2,9,2),NF
   23 GOTO (2,11,2,2,2,2,2,9,2),NF
    9 SZ=CMPLX(S,-PI)
      GOTO 2
   10 SZ=CMPLX(S,0.0)
      GOTO 2
   11 SZ=CMPLX(S,PI)
    2 CLOGOK=SZ
C*UL  80  RETURN
      RETURN
      END
