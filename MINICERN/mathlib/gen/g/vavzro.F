*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE VAVZRO(A,B,X,RKA,BE2,LU)
C
C     RZERO SEARCHES FOR THE ROOT OF THE EQUATION FCN=0 IN THE INTERVAL
C
      COMMON /FORFCN/ SS,LFCN
      DATA E,EPSI,MAXFUN/1E-9,1E-5,100/
      MC=0
C
      XA=MIN(A,B)
      XB=MAX(A,B)
      FA=VAVFCN(A,RKA,BE2)
      MC=MC+1
      FB=VAVFCN(B,RKA,BE2)
      IF(FA*FB .GT. 0.0) GO TO 16
      MC=MC+1
C
    4 X=0.5*(XA+XB)
      R=X-XA
      EE=ABS(X)+E
      IF(R .LE. EE*EPSI) GO TO 18
      F1=FA
      X1=XA
      F2=FB
      X2=XB
    1 MC=MC+1
      G=VAVFCN(X,RKA,BE2)
      IF(MC .GT. MAXFUN) GO TO 17
      FX=G
C
      IF(FX*FA .GT. 0.0) GO TO 2
      FB=FX
      XB=X
      GO TO 3
    2 XA=X
      FA=FX
C
C     PARABOLA ITERATION
C
    3 IF((X1-X2)*(X2-X)*(X1-X) .EQ. 0.0) GO TO 4
      F3=FX
      X3=X
      U1=(F1-F2)/(X1-X2)
      U2=(F2-FX)/(X2-X)
      CA=U1-U2
      CB=(X1+X2)*U2-(X2+X)*U1
      CC=(X1-X)*F1-X1*(CA*X1+CB)
      IF(CA .EQ. 0.0) GO TO 8
      U3=0.5*CB/CA
      U4=U3**2-CC/CA
      IF(U4 .LT. 0.0) GO TO 4
      U5=SQRT(U4)
      IF(X .GE. -U3) GO TO 10
      X=-U3-U5
      GO TO 9
   10 X=-U3+U5
      GO TO 9
    8 X=-CC/CB
    9 IF(X .LT. XA) GO TO 4
      IF(X .GT. XB) GO TO 4
C
C     TEST FOR OUTPUT
C
      R=ABS(X-X3)
      R1=ABS(X-X2)
      IF(R .GT. R1) R=R1
      EE=ABS(X)+E
      IF(R/EE .GT. EPSI) GO TO 5
      MC=MC+1
      G=VAVFCN(X,RKA,BE2)
      IF(MC .GT. MAXFUN) GO TO 17
      FX=G
      IF(FX .EQ. 0.0) GO TO 18
      IF(FX*FA .LT. 0.0) GO TO 7
      XX=X+EPSI*EE
      IF(XX .GE. XB) GO TO 18
      MC=MC+1
      G=VAVFCN(X,RKA,BE2)
      IF(MC .GT. MAXFUN) GO TO 17
      FF=VAVFCN(XX,RKA,BE2)
      FA=FF
      XA=XX
      GO TO 6
    7 XX=X-EPSI*EE
      IF(XX .LE. XA) GO TO 18
      MC=MC+1
      FX=G
      IF(MC .GT. MAXFUN) GO TO 17
      FF=VAVFCN(XX,RKA,BE2)
      FB=FF
      XB=XX
    6 IF(FX*FF .GT. 0.0) GO TO 14
   18 R=EPSI*EE
      FF=VAVFCN(X,RKA,BE2)
      RETURN
   14 F1=F2
      X1=X2
      F2=FX
      X2=X
      X=XX
      FX=FF
      GO TO 3
C
    5 F1=F2
      X1=X2
      F2=F3
      X2=X3
      GO TO 1
C
   16 WRITE(LU,301)
      RETURN
C
   17 WRITE(LU,300) X,G,LFCN
  301 FORMAT(/10X,' RZERO   FCN(A) AND FCN(B) HAVE SAME SIGN'/)
  300 FORMAT(/10X,' RZERO   NUMBER OF ITERATIONS EXCEEDED'/
     1  10X,' X=',E15.5,2X,' FCN(X)=',E15.5,2X,' LFCN=',I2/)
      RETURN
C
      END
