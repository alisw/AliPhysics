*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SGZ (M,S)
C-----SUBROUTINE SGZ SETS UP S(OBJECTIVE FUNCTION),G(GRADIENT OF S) AND
C-----Z(APPROXIMATE COVARIANCE MATRIX)
#include "d510si.inc"
#include "d510ui.inc"
      K=NED(1)
      K2=1
      DO 12 L1=1,K
      K1=K2
      NX=NED(2)-2
      IF (INDFLG(3)) 1,2,1
 1    NX=NED(2)
      K1=K1-2
 2    CONTINUE
      DO 3 I=1,NX
      KI=K1+1+I
      X(I)=EXDA(KI)
 3    CONTINUE
      CALL ARITHM (Y)
      IF (INDFLG(3)) 4,6,4
 4    IF (Y) 13,13,5
C-----MAXIMUM LIKLEHOOD
 5    S=S-LOG(Y)
      Y=-Y
      SIG=Y
      GO TO 7
C-----CHI SQUARED
 6    SIG=EXDA(K2+1)
      Y=Y-EXDA(K1)
      S=S+((Y/SIG)**2)/2.
 7    CONTINUE
      N=0
      DO 9 J=1,M
      IF (PL0(J)) 9,9,8
 8    N=N+1
      DF(N)=DF(J)/SIG
      G(J)=G(J)+DF(N)*(Y/SIG)
 9    CONTINUE
      L=1
      IF (N.LT.1) GO TO 11
      DO 10 I=1,N
      DO 10 J=1,I
      Z(L)=Z(L)+DF(I)*DF(J)
      L=L+1
 10   CONTINUE
 11   K2=K2+NED(2)
 12   CONTINUE
      RETURN
C-----  -VE OR ZERO Y IN MAXIMUM LIKLEHOOD
 13   INDFLG(1)=1
      S=1.E10
      RETURN
      END
