*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION VVIDEN(X)

      LOGICAL LST

      COMMON /G116C1/ H(7),T0,T1,T,OMEGA,A(155),B(155),X0

      PARAMETER (PI = 3.14159 265, Z1 = 1, HF = Z1/2)

      LST=.FALSE.
      GO TO 9

      ENTRY VVIDIS(X)

      LST=.TRUE.
    9 N=X0
      IF(X .LT. T0) THEN
       F=0
      ELSEIF(X .LE. T1) THEN
       Y=X-T0
       U=OMEGA*Y-PI
       COF=2*COS(U)
       A1=0
       A0=A(1)
       DO 1 K = 2,N+1
       A2=A1
       A1=A0
    1  A0=A(K)+COF*A1-A2
       B1=0
       B0=B(1)
       DO 2 K = 2,N
       B2=B1
       B1=B0
    2  B0=B(K)+COF*B1-B2
       F=HF*(A0-A2)+B0*SIN(U)
       IF(LST) F=F+Y/T
      ELSE
       F=0
       IF(LST) F=1
      ENDIF
      VVIDEN=F
      RETURN
      END
