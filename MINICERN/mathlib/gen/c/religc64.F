*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:11  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DELIGC(AKP,P,A,B)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELIGC/DELIGC')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RELIGC(AKP,P,A,B)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELIGC')
#endif
C
C     Translation of Algol procedure cel(kc,p,a,b) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions III.,  Numer. Math. 13 (1969) 305-315
C
      PARAMETER (ID = 16)
      PARAMETER (PI = 3.14159 26535 89793 24D0, PIH = PI/2)
      PARAMETER (Z10 = 10)
      PARAMETER (CA = Z10**(-ID/2))

      IF(AKP .EQ. 0) THEN
       H=0
       CALL MTLPRT(NAME,'C347.4','AKP = 0')
      ELSE
       PP=P
       AA=A
       BB=B
       YKP=ABS(AKP)
       E=YKP
       XM=1
       IF(PP .GT. 0) THEN
        PP=SQRT(PP)
        BB=BB/PP
       ELSE
        F=YKP**2
        Q=1-F
        G=1-PP
        F=F-PP
        Q=(BB-AA*PP)*Q
        PP=SQRT(F/G)
        AA=(AA-BB)/G
        BB=-Q/(G**2*PP)+AA*PP
       ENDIF
    1  F=AA
       AA=BB/PP+AA
       G=E/PP
       BB=2*(F*G+BB)
       PP=G+PP
       G=XM
       XM=YKP+XM
       IF(ABS(G-YKP) .GT. G*CA) THEN
        YKP=2*SQRT(E)
        E=YKP*XM
        GO TO 1
       ENDIF
       H=PIH*(AA*XM+BB)/(XM*(XM+PP))
      ENDIF
#if defined(CERNLIB_DOUBLE)
      DELIGC=H
#endif
#if !defined(CERNLIB_DOUBLE)
      RELIGC=H
#endif
      RETURN
      END
