*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_SINGLE)
      FUNCTION C309R4(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     1              FPMIN,FPMAX,LPR)
C
C ***   Evaluate CF1 = F'(ZL,ETA,X)/F(ZL,ETA,X)    (REAL)
C
      IMPLICIT REAL(A-H,O-Z)
      LOGICAL LPR,ETANE0

      FCL=1
      XI=1/X
      PK=ZL+1
      PX=PK+LIMIT
      EK=ETA/PK
      F=EK+PK*XI
      IF(ABS(F) .LT. FPMIN) F=FPMIN
      D=0
      C=F
      SMALL=SQRT(FPMIN)
      RK2=1+EK*EK
C
C ***   begin CF1 loop on PK = k = lambda + 1
C
   10 PK1=PK+1
      TPK1=PK+PK1
      IF(ETANE0) THEN
       EK=ETA/PK
       RK2=1+EK*EK
       TK=TPK1*(XI+EK/PK1)
      ELSE
       TK=TPK1*XI
      END IF
      C=TK-RK2/C
      D=TK-RK2*D
      IF(ABS(C) .LT. FPMIN) C=FPMIN
      IF(ABS(D) .LT. FPMIN) D=FPMIN
      D=1/D
      DF=D*C
      F=F*DF
      FCL=FCL*D*TPK1*XI
      IF(ABS(FCL) .LT. SMALL) FCL=FCL/SMALL
      IF(ABS(FCL) .GT. FPMAX) FCL=FCL*FPMIN
      PK=PK1
      IF(PK .LE. PX) THEN
       IF(ABS(DF-1) .GE. EPS) GO TO 10
       NFP=PK-ZL-1
       ERR=EPS*SQRT(REAL(NFP))
       C309R4=F
      ELSE
       IF(LPR) WRITE (6,1000) LIMIT,ABS(X)
       ERR=2
      END IF
      RETURN
 1000 FORMAT(1X,'***** CERN C309 CCLBES ... CF1 (REAL) HAS FAILED ',
     1'TO CONVERGE AFTER',I10,' ITERATIONS AS ABS(X) =',F15.0)
      END
#endif
