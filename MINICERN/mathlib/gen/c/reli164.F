*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:10  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DELI1(X,AKP)
C
#include "gen/imp64.inc"
C
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RELI1(X,AKP)
C
#endif
C
C     Translation of Algol procedure el1(x,kc) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions,  Numer. Math. 7 (1965) 78-90
C
      PARAMETER (ID = 16)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (Z10 = 10)
      PARAMETER (CA = Z10**(-ID/2), CB = Z10**(-(ID+2)))

      IF(X .EQ. 0) THEN
       H=0
      ELSEIF(AKP .EQ. 0) THEN
#if defined(CERNLIB_DOUBLE)
       H=DASINH(X)
#endif
#if !defined(CERNLIB_DOUBLE)
       H=ASINH(X)
#endif
      ELSE
       Y=ABS(1/X)
       YKP=ABS(AKP)
       XM=1
       L=0
    1  E=XM*YKP
       G=XM
       XM=YKP+XM
       Y=-E/Y+Y
       IF(Y .EQ. 0) Y=SQRT(E)*CB
       IF(ABS(G-YKP) .GT. CA*G) THEN
        YKP=2*SQRT(E)
        L=2*L
        IF(Y .LT. 0) L=L+1
        GO TO 1
       ENDIF
       IF(Y .LT. 0) L=L+1
       E=(ATAN(XM/Y)+PI*L)/XM
       IF(X .LT. 0) E=-E
       H=E
      ENDIF
#if defined(CERNLIB_DOUBLE)
      DELI1=H
#endif
#if !defined(CERNLIB_DOUBLE)
      RELI1=H
#endif
      RETURN
      END
