*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:47  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE COEDIS(RKA,BE2,I,J)
C
C     COEDIS COMPUTES THE ENDPOINTS T- AND T+ OF THE SUPPORT OF
C     DISVAV(X,0).IT ALSO COMPUTES THE FOURIER COEFFICIENTS OF DISVAV(X,
C
      COMMON /VAVILI/ T0,T1,T,OMEGA
      COMMON /VAVILO/ A(155),B(155),N
      COMMON /FORFCN/ SS,LFCN
      DIMENSION XP(8),XQ(6)
      DATA E,PI,RG /5E-4, 3.1415926535898, 0.5772156649015/
      DATA  XP
     + /9.29, 2.47, 0.89, 0.36, 0.15, 0.07, 0.03, 0.02/
      DATA  XQ
     + /0.012, 0.03, 0.08, 0.26, 0.87, 3.83/
C
      LU=ABS(J)
      IF (RKA .LT. 0.01 .OR. RKA .GT. 10.0) GOTO 6
      IF (BE2 .LT. 0.0  .OR. BE2 .GT.  1.0) GOTO 8
C
      Z=1.-BE2*(1.-RG)-LOG(E)/RKA
      T0=(LOG(E)/RKA-(1.+BE2*RG)-Z*LOG(RKA)-(Z+BE2)*(LOG(Z)
     1  +EXPINT(Z))+EXP(-Z))/Z
      DO 1 L = 1,8
      IF(RKA .GE. XP(L)) GO TO 11
    1 CONTINUE
      L=9
   11 P=-L-0.5
      DO 2 L = 1,6
      IF(RKA .LE. XQ(L)) GO TO 22
    2 CONTINUE
      L=7
   22 Q=L-7.5
      LFCN=3
      CALL VAVZRO(P,Q,X,RKA,BE2,LU)
      T1=(LOG(E)/RKA-(1.+BE2*RG))/X-LOG(RKA)-(1.+BE2/X)*(LOG(ABS(X))
     1  +EXPINT(X))+EXP(-X)/X
C
      IF(J .GT. 0) WRITE(J,10) T0,T1
      T=T1-T0
      OMEGA=2.0*PI/T
      LFCN=1
      CALL VAVZRO(5.,155.,X,RKA,BE2,LU)
      N=X+1.
C
      D=EXP(RKA*(1.+BE2*(RG-LOG(RKA))))
      A(N)=0.
      IF(I .EQ. 0) A(N)=OMEGA/PI
      N1=N-1
      Q=-1.
C
      DO 3 K = 1,N1
      L=N-K
      X=OMEGA*K
      X1=X/RKA
      C1=LOG(X)-COSINT(X1)
      C2=SININT(X1)
      C3=SIN(X1)
      C4=COS(X1)
      F1=RKA*(BE2*C1-C4)-X*C2
      F2=X*C1+RKA*(C3+BE2*C2)+T0*X
      D1=Q*D*EXP(F1)/PI
      HS=D1*SIN(F2)
      HC=D1*COS(F2)
      IF(I .EQ. 0) GO TO 4
      A(L)=HS/K
      B(L)=HC/K
      A(N)=A(N)-2.0*Q*A(L)
      GO TO 3
    4 A(L)=HC*OMEGA
      B(L)=-HS*OMEGA
    3 Q=-Q
      RETURN
C
    6 WRITE(LU,7) RKA
      RETURN
C
    8 WRITE(LU,9) BE2
    7 FORMAT(/10X,'KAPPA ='  ,E10.3,' - OUT OF RANGE')
    9 FORMAT(/10X,'BETA**2 =',E10.3,' - OUT OF RANGE')
   10 FORMAT(10X,'T- =',F8.3,10X,'T+ =',F8.3)
      RETURN
C
      END
