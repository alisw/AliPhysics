*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R1(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH,
     1             LPR,ACCUR,DELL)
C
C                                    (omega)        (omega)
C *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
C                                    ZL             ZL
C     where PM = omega.i
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL LPR
      DOUBLE PRECISION EPS,ERR,ACC8,ACCH,ACCUR,TA,RK
      DOUBLE PRECISION ABSC,HALF

      PARAMETER(HALF = 1D0/2D0)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))

      TA=LIMIT+LIMIT
      ETAP=ETA*PM
      XI=1/X
      WI=ETAP+ETAP
      RK=0
      PQ=(1-ETA*XI)*PM
      AA=-(ETA*ETA+ZL*ZL+ZL)+ETAP
      BB=2*(X-ETA+PM)
      RL=XI*PM
      IF(ABSC(BB) .LT. ACCH) THEN
       RL=RL*AA/(AA+RK+WI)
       PQ=PQ+RL*(BB+PM+PM)
       AA=AA+2*(RK+1+WI)
       BB=BB+4*PM
       RK=RK+4
      END IF
      DD=1/BB
      DL=AA*DD*RL
   10 PQ=PQ+DL
      RK=RK+2
      AA=AA+RK+WI
      BB=BB+PM+PM
      DD=1/(AA*DD+BB)
      DL=DL*(BB*DD-1)
      ERR=ABSC(DL)/ABSC(PQ)
      IF(ERR .GE. MAX(EPS,ACC8*RK*HALF) .AND. RK .LE. TA) GO TO 10
C
      NPQ=HALF*RK
      C309R1=PQ+DL
      IF(LPR .AND. NPQ .GE. LIMIT-1 .AND. ERR .GT. ACCUR)
     1 WRITE(6,1000) INT(DIMAG(PM)),NPQ,ERR,ZL+DELL
      RETURN
 1000 FORMAT(1X,'***** CERN C309 WCLBES ... ',
     2 'CF2(',I2,') NOT CONVERGED FULLY IN ',I7,' ITERATIONS'/1X,27X,
     3 'ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL = ',2F8.3)
      END
#endif
