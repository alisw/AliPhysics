*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE H101S1(A,B,C,Z,M,N,IDA,IP,KP,LW,IDW,EPS)

#include "gen/imp64.inc"

      DIMENSION A(IDA,*),B(*),C(*),LW(IDW,*)

      PARAMETER (R10 = 10, IE0 = 15)

C     Exchanges a basic with a non-basic variable and transforms
C     the whole tableau

      A(IP,KP)=1/A(IP,KP)
      DO 1 I = 1,M
      IF(I .NE. IP) A(I,KP)=A(I,KP)*A(IP,KP)
    1 CONTINUE
      C(KP)=C(KP)*A(IP,KP)
      DO 3 I = 1,M
      IF(I .NE. IP) THEN
       DO 2 K = 1,N
       IF(K .NE. KP) A(I,K)=A(I,K)-A(IP,K)*A(I,KP)
    2  CONTINUE
       B(I)=B(I)-B(IP)*A(I,KP)
      ENDIF
    3 CONTINUE
      DO 4 K = 1,N
      IF(K .NE. KP) THEN
       C(K)=C(K)-A(IP,K)*C(KP)
       A(IP,K)=-A(IP,K)*A(IP,KP)
      ENDIF
    4 CONTINUE
      Z=Z-B(IP)*C(KP)
      B(IP)=-B(IP)*A(IP,KP)

      EPS=0
      DO 5 I = 1,M
      DO 5 K = 1,N
    5 EPS=EPS+ABS(A(I,K))
      EPSL=LOG10(2*EPS/(M*N))
      IEXP=INT(EPSL)-IE0
      IF(EPSL .LT. 0) IEXP=IEXP-1
      EPS=R10**IEXP

      DO 10 I = 1,M
      IF(ABS(B(I)) .LT. EPS) B(I)=0
      DO 10 K = 1,N
      IF(ABS(A(I,K)) .LT. EPS) A(I,K)=0
   10 CONTINUE
      DO 25 K = 1,N
      IF(ABS(C(K)) .LT. EPS) C(K)=0
   25 CONTINUE
      IF(ABS(Z) .LT. EPS) Z=0
      IR=LW(IP,4)
      LW(IP,4)=LW(KP,5)
      LW(KP,5)=IR
      RETURN
      END
