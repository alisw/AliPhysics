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
      FUNCTION DELI1C(AKP)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME1,NAME2,NAME3
      PARAMETER(NAME1='RELI1C/DELI1C')
      PARAMETER(NAME2='RELI2C/DELI2C')
      PARAMETER(NAME3='RELI3C/DELI3C')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RELI1C(AKP)
C
      CHARACTER*(*) NAME1,NAME2,NAME3
      PARAMETER(NAME1='RELI1C')
      PARAMETER(NAME2='RELI2C')
      PARAMETER(NAME3='RELI3C')
#endif
C
C     Translation of Algol procedures cel1(kc), cel2(kc,a,b) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions,  Numer. Math. 7 (1965) 78-90
C     and of Algol procedure cel3(kc,m,p) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions II.,  Numer. Math. 7 (1965) 353-354
C
      PARAMETER (ID = 16)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PIH = PI/2, PIQ = PI/4)
      PARAMETER (Z1 = 1, Z10 = 10, HF = Z1/2)
      PARAMETER (CA = Z10**(-ID/2))

      IF(AKP .EQ. 0) THEN
       H=0
       CALL MTLPRT(NAME1,'C347.1','AKP = 0')
      ELSE
       YKP=ABS(AKP)
       XM=1
    1  G=XM
       XM=YKP+XM
       IF(ABS(G-YKP) .GT. CA*G) THEN
        YKP=SQRT(G*YKP)
        XM=HF*XM
        GO TO 1
       ENDIF
       H=PI/XM
      ENDIF
      GO TO 9

#if defined(CERNLIB_DOUBLE)
      ENTRY DELI2C(AKP,A,B)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY RELI2C(AKP,A,B)
#endif

      IF(AKP .EQ. 0) THEN
       IF(B .EQ. 0) THEN
        H=A
       ELSE
        H=0
        CALL MTLPRT(NAME2,'C347.2','AKP = 0, B NE 0')
       ENDIF
      ELSE
       AA=A
       BB=B
       XM=1
       C=AA
       AA=BB+AA
       YKP=ABS(AKP)
    2  BB=2*(C*YKP+BB)
       C=AA
       XM0=XM
       XM=YKP+XM
       AA=BB/XM+AA
       IF(ABS(XM0-YKP) .GT. CA*XM0) THEN
        YKP=2*SQRT(YKP*XM0)
        GO TO 2
       ENDIF
       H=PIQ*AA/XM
      ENDIF
      GO TO 9

#if defined(CERNLIB_DOUBLE)
      ENTRY DELI3C(AKP,AK2,P)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY RELI3C(AKP,AK2,P)
#endif

      IF(AKP*P .EQ. 0) THEN
       H=0
       CALL MTLPRT(NAME3,'C347.3','AKP * P = 0')
      ELSE
       YKP=ABS(AKP)
       E=YKP
       AM0=1
       PP=P
       IF(PP .GT. 0) THEN
        C=1
        PP=SQRT(PP)
        D=1/PP
       ELSE
        G=1-PP
        F=YKP**2-PP
        PP=SQRT(F/G)
        D=-AK2/(G*PP)
        C=0
       ENDIF
    3  F=C
       C=D/PP+C
       G=E/PP
       D=2*(F*G+D)
       PP=G+PP
       G=AM0
       AM0=YKP+AM0
       IF(ABS(G-YKP) .GT. CA*G) THEN
        YKP=2*SQRT(E)
        E=YKP*AM0
        GO TO 3
       ENDIF
       H=PIH*(C*AM0+D)/(AM0*(AM0+PP))
      ENDIF
#if defined(CERNLIB_DOUBLE)
    9 DELI1C=H
#endif
#if !defined(CERNLIB_DOUBLE)
    9 RELI1C=H
#endif
      RETURN
      END
