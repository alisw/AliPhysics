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
      FUNCTION DELI2(X,AKP,A,B,MODE)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELI2/DELI2')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RELI2(X,AKP,A,B,MODE)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RELI2')
#endif
C
C     Translation of Algol procedure el2(x,kc,a,b) in
C      R. BULIRSCH Numerical Calculation of Elliptic Integrals and
C      Elliptic Functions,  Numer. Math. 7 (1965) 78-90,
C      extended for (k sub c)**2 < 0  (MODE = -1)
C
      PARAMETER (ID = 16)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (Z1 = 1, Z10 = 10, HF = Z1/2)
      PARAMETER (CA = Z10**(-ID/2), CB = Z10**(-(ID+2)))

      CHARACTER*80 ERRTXT

      IF(X .EQ. 0 .OR. AKP .EQ. 0) THEN
#if defined(CERNLIB_DOUBLE)
       H=B*DASINH(X)+(A-B)*X/SQRT(1+X**2)
#endif
#if !defined(CERNLIB_DOUBLE)
       H=B*ASINH(X)+(A-B)*X/SQRT(1+X**2)
#endif
      ELSE
       IF(MODE .EQ. 1) THEN
        XX=X
        BKP=AKP
        AA=A
        BB=B
       ELSEIF(MODE .EQ. -1) THEN
        W=1-(AKP*X)**2
        IF(W .GT. 0) THEN
         W1=1/SQRT(1+AKP**2)
         XX=X/(W1*SQRT(W))
         BKP=AKP*W1
         AA=A*W1
         BB=(B+A*AKP**2)*W1**3
        ELSE
         H=0
         WRITE(ERRTXT,101) X,AKP
         CALL MTLPRT(NAME,'C346.1',ERRTXT)
         GO TO 9
        ENDIF
       ELSE
        H=0
        WRITE(ERRTXT,102) MODE
        CALL MTLPRT(NAME,'C346.2',ERRTXT)
        GO TO 9
       ENDIF
       C=XX**2
       D=1+C
       PP=SQRT((1+BKP**2*C)/D)
       D=XX/D
       C=D/(2*PP)
       Z=AA-BB
       XI=AA
       AA=HF*(AA+BB)
       Y=ABS(1/XX)
       F=0
       L=0
       XM=1
       YKP=ABS(BKP)
    1  BB=XI*YKP+BB
       E=XM*YKP
       G=E/PP
       D=F*G+D
       F=C
       XI=AA
       PP=G+PP
       C=HF*(D/PP+C)
       G=XM
       XM=YKP+XM
       AA=HF*(BB/XM+AA)
       Y=Y-E/Y
       IF(Y .EQ. 0) Y=SQRT(E)*CB
       IF(ABS(G-YKP) .GT. CA*G) THEN
        YKP=2*SQRT(E)
        L=2*L
        IF(Y .LT. 0) L=L+1
        GO TO 1
       ENDIF
       IF(Y .LT. 0) L=L+1
       E=(ATAN(XM/Y)+PI*L)*AA/XM
       IF(XX .LT. 0) E=-E
       H=E+C*Z
      ENDIF
#if defined(CERNLIB_DOUBLE)
    9 DELI2=H
#endif
#if !defined(CERNLIB_DOUBLE)
    9 RELI2=H
#endif
      RETURN
  101 FORMAT('X = ',1P,D15.8,' AKP = ',D15.8,' ILLEGAL',
     1       '[(AKP * X)**2 >= 1]')
  102 FORMAT('MODE = ',I5,' ILLEGAL')
      END
