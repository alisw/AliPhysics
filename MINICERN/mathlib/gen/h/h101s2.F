*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE H101S2(A,B,C,M,N,IDA,IP,KP,KK,LW,IDW,W,X,EPS,IRC)

#include "gen/imp64.inc"

      DIMENSION A(IDA,*),B(*),C(*),X(*),W(*),LW(IDW,*)

C     Finds the pivot-column, by taking degeneracy into account.
C     Columns with C(K) >= 0 and A(IP,K) > 0: LW(.,3) = 1

      DO 1 K = 1,N
    1 LW(K,3)=0
      DO 5 K = 1,N
      IF(LW(K,2) .EQ. IRC .AND. K .NE. KK) THEN
       IF(ABS(A(IP,K)) .LT. EPS) A(IP,K)=0
       IF(ABS(C(K)) .LT. EPS) C(K)=0
       IF(A(IP,K) .GT. 0 .AND. C(K) .GE. 0) THEN
        LW(K,3)=1
        Q=C(K)/A(IP,K)
       ENDIF
      ENDIF
    5 CONTINUE
      KP=0
      DO 6 K = 1,N
    6 KP=KP+LW(K,3)
      IF(KP .EQ. 0) RETURN

C     KP = 0 : No column found.
C     Only columns with minimum quotient: LW(.,3) = 1.

      DMIN=Q
      DO 11 J = 1,2
      DO 11 K = 1,N
      IF(LW(K,3) .NE. 0) THEN
       Q=C(K)/A(IP,K)
       IF(Q .LE. DMIN) THEN
        DMIN=Q
        KP=K
        GO TO 11
       ENDIF
       LW(K,3)=0
      ENDIF
   11 CONTINUE

   12 IND=0
      DO 15 K = 1,N
   15 IND=IND+LW(K,3)
      IF(IND .EQ. 1) RETURN

C     Two possible columns are picked out

      KP1=0
      DO 20 K = 1,N
      IF(LW(K,3) .NE. 0) THEN
       IF(KP1 .EQ. 0) KP1=K
       KP2=K
      ENDIF
   20 CONTINUE

C     Choose pivot column from two columns with equal quotient.

      DO 25 J = 1,M+N
      X(J)=0
   25 W(J)=0
      X(LW(KP1,5))=1
      W(LW(KP2,5))=1
      DO 30 I = 1,M
      IF(LW(I,1) .EQ. IRC) THEN
       IF(ABS(A(I,KP1)) .LT. EPS) A(I,KP1)=0
       IF(ABS(A(I,KP2)) .LT. EPS) A(I,KP2)=0
       X(LW(I,4))=A(I,KP1)
       W(LW(I,4))=A(I,KP2)
      ENDIF
   30 CONTINUE
      DO 35 J = 1,M+N
      IF(J .EQ. LW(IP,4)) THEN
       X(J)=1
       W(J)=1
      ELSE
       X(J)=X(J)/A(IP,KP1)
       W(J)=W(J)/A(IP,KP2)
      ENDIF
   35 CONTINUE
      DO 50 J = 1,M+N
      IF(X(J) .LT. W(J)) THEN
       LW(KP2,3)=0
       KP=KP1
       GO TO 12
      ELSEIF(X(J) .GT. W(J)) THEN
       LW(KP1,3)=0
       KP=KP2
       GO TO 12
      ENDIF
   50 CONTINUE
      RETURN
      END
