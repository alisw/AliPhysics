*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION  CHSUM(MODE,C,N,X)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DCHSUM(MODE,C,N,X)
#endif
#include "gen/imp64.inc"
      DIMENSION C(0:N)

      IF(MODE .EQ. 1) THEN
       H=X
       F=H
       V=1
      ELSE IF(MODE .EQ. 2) THEN
       H=2*X**2-1
       F=H
       V=1
      ELSE IF(MODE .EQ. 3) THEN
       H=2*X**2-1
       F=1
       V=X
      ELSE IF(MODE .EQ. 4) THEN
       H=2*X-1
       F=H
       V=1
      END IF

      ALFA=H+H
      B1=0
      B2=0
      DO 1 I = N,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
#if !defined(CERNLIB_DOUBLE)
       CHSUM=V*(B0-F*B2)
#endif
#if defined(CERNLIB_DOUBLE)
      DCHSUM=V*(B0-F*B2)
#endif
      RETURN
      END
