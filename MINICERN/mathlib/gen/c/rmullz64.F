*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DMULLZ(A,N,MAXITR,Z0)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RMULLZ(A,N,MAXITR,Z0)
#endif
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +  Z0,Z,DX,X,X3,Y1,Y2,Y,TE(7)
      LOGICAL LSW
      DIMENSION A(0:*),Z0(*)
      CHARACTER*(*) NAME
#if defined(CERNLIB_DOUBLE)
      PARAMETER( NAME='DMULLZ')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER( NAME='RMULLZ')
#endif
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (ETA1 = 1D-14, ETA2 = 6D-8, BIG = 1D20)
      PARAMETER (C1 = 0.9D0, C2 = 0.002D0, C3 = 0.1D0)

#include "gen/gcmpfun.inc"
      SUMABS(Z)=ABS(GREAL(Z))+ABS(GIMAG(Z))

      IF(N .LE. 0) RETURN
      IF(A(0) .EQ. 0) THEN
       CALL MTLPRT(NAME,'C202.1','A(0) = 0')
       RETURN
      ENDIF

      N1=N
    2 IF(N1 .EQ. 1) THEN
       Z0(1)=-A(1)/A(0)
       RETURN
      ENDIF
      IF(A(N1) .EQ. 0) THEN
       Z0(N1)=0
       N1=N1-1
       GO TO 2
      ENDIF

      SCALE=ABS(A(N1)/A(0))**(Z1/N1)
      B=A(0)
      DO 6 I = 1,N1
      B=B*SCALE
    6 Z0(I)=A(I)/B
      IF(N1 .EQ. 2) GO TO 1

   10 LSW=.TRUE.
      Y1=Z0(1)+1
      Y2=Z0(1)-1
      DO 11 I = 2,N1
      Y1=Z0(I)+Y1
   11 Y2=Z0(I)-Y2
      Y=Z0(N1)
      X=0
      DX=1

      TE(1)=-2
   12 TE(2)=Y2/Y
      TE(3)=(Y1-Y2)/(Y*TE(1))
      DO 17 ITER = 1,MAXITR
      TE(4)=TE(2)-1
      TE(5)=(TE(4)-TE(3))/(TE(1)+1)
      TE(6)=HF*(TE(5)+TE(4))
      TE(7)=SQRT(TE(6)**2+TE(5))
      TE(1)=TE(6)+TE(7)
      TE(7)=TE(6)-TE(7)
      IF(ABS(TE(1)) .LE. ABS(TE(7))) THEN
       IF(TE(7) .EQ. 0) TE(7)=C1
       TE(1)=TE(7)
      ENDIF
      DX=DX/TE(1)
      X=DX+X
      EPSI=SUMABS(X)*ETA1
      IF(SUMABS(DX) .LT. EPSI .AND. SUMABS(Y) .LT. C2) GO TO 18
      Y2=Y
      Y=X+Z0(1)
      DO 21 I = 2,N1
   21 Y=Y*X+Z0(I)
      IF(Y .EQ. 0) GO TO 18

   15 IF(SUMABS(Y) .GE. 100*SUMABS(Y2) .AND. SUMABS(DX) .GE. EPSI) THEN
       TE(1)=2*TE(1)
       DX=HF*DX
       X=X-DX
       Y=X+Z0(1)
       DO 22 I = 2,N1
   22  Y=Y*X+Z0(I)
       IF(Y .EQ. 0) GO TO 18
       GOTO 15
      ENDIF

      TE(2)=Y2/Y
      TE(3)=(TE(2)/TE(1))*TE(4)
   17 CONTINUE

      CN=ABS(Z0(N1))
      IF(ABS(CN-1) .LT. C3) THEN
       DO 40 I = 1,N1
   40  Z0(I)=BIG
       CALL MTLPRT(NAME,'C202.2','TOO MANY ITERATIONS')
       RETURN
      ENDIF
      S=CN**(-Z1/N1)
      SCALE=SCALE/S
      B=1
      DO 30 I = 1,N1
      B=B*S
   30 Z0(I)=Z0(I)*B
      GO TO 10

   20 IF(ABS(GIMAG(X)) .LT. ABS(GREAL(X))*ETA2) GO TO 10
      LSW=.FALSE.
      X3=GCONJG(X)
      DX=GCONJG(DX)
      TE(1)=GCONJG(TE(1))
      X=X3-DX
      Y=X+Z0(1)
      DO 51 I = 2,N1
   51 Y=Y*X+Z0(I)
      IF(Y .EQ. 0) GO TO 18
      Y2=Y
      X=X-DX*TE(1)
      Y=X+Z0(1)
      DO 52 I = 2,N1
   52 Y=Y*X+Z0(I)
      IF(Y .EQ. 0) GO TO 18
      Y1=Y
      X=X3
      Y=X+Z0(1)
      DO 53 I = 2,N1
   53 Y=Y*X+Z0(I)
      IF(Y .NE. 0) GO TO 12

   18 Z0(N1)=X*SCALE
      N1=N1-1
      Z0(1)=X+Z0(1)
      DO 19 I = 2,N1
   19 Z0(I)=Z0(I-1)*X+Z0(I)
      IF(N1 .GT. 2) THEN
       IF(LSW) GO TO 20
       GO TO 10
      ENDIF

    1 Z0(2)=(SQRT((HF*Z0(1))**2-Z0(2))-HF*Z0(1))*SCALE
      Z0(1)=-Z0(1)*SCALE-Z0(2)
      RETURN
      END
