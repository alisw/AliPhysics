*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION STUDIN(Q,N)
C    THIS ROUTINE COMPUTES THE INVERSE OF THE DISTRIBUTION
C    FUNCTION OF THE STUDENT DISTRIBUTION WITH N DEGREES OF FREEDOM.
C    Q=PROBABILITY AT WHICH THE FUNCTION IS INVERTED,
C    N=INTEGER GREATER THAN ZERO.
C     THE ROUTINE WAS WRITTEN BY G.W.HILL IN ALGOL
C     C.A.C.M. ALGORITHM 396
C    AT LEAST 6 SIGNIFICANT FIGURES ARE CORRECT.
      DATA HP/1.5707963268/
      IF(N. LT. 1) GO TO 10
      RL=1.
      IF(Q .GE. 0.5) GO TO 1
      RL=-1.
      P=2.*Q
      GO TO 2
    1 P=2.*(1.-Q)
    2 IF(P .LE. 0. .OR. P .GT. 1.) GO TO 20
      IF(N .GT. 1) GO TO 3
      PP=COS(HP*P)
      STUDIN=PP/SQRT(1.-PP*PP)*RL
      RETURN
    3 IF(N .GT. 2) GO TO 4
      STUDIN=SQRT(2./(P*(2.-P))-2.)*RL
      RETURN
    4 RN=N
      A=1./(RN-0.5)
      B=48./(A*A)
      C=((20700.*A/B-98.)*A-16.)*A+96.36
      D=((94.5/(B+C)-3.)/B+1.)*SQRT(A*HP)*RN
      X=D*P
      Y=X**(2./RN)
      IF(Y .LE. 0.05+A) GO TO 5
      PP=0.5*P
      X=GAUSIN(PP)
      Y=X*X
      IF(N .GE. 5) C=C+0.3*(RN-4.5)*(X+0.6)
      C=(((0.05*D*X-5.)*X-7.)*X-2.)*X+B+C
      Y=(((((0.4*Y+6.3)*Y+36.)*Y+94.5)/C-Y-3.)/B+1.)*X
      Y=A*Y*Y
      IF(Y .LE. 0.002) Y=0.5*Y*Y+Y
      IF(Y .GT. 0.002) Y=EXP(Y)-1.
      GO TO 6
    5 Y=((1./(((RN+6.)/(RN*Y)-0.089*D-0.822)*(RN+2.)*3.)+0.5/(RN+4.))*Y-
     11.)*(RN+1.)/(RN+2.)+1./Y
    6 STUDIN=SQRT(RN*Y)*RL
      RETURN
   10 WRITE(6,7) N
   20 WRITE(6,8) Q
      STOP
    7 FORMAT(/10X,'DEGREE OF FREEDOM N=',I5,' IN STUDIN ILLEGAL'/)
    8 FORMAT(/10X,'ARGUMENT  Q=',E15.5,' IN STUDIN ILLEGAL'/)
      END
