*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE NMDCHL (N,NADIM,AHESS,EPSMCH,Z,P)
      INTEGER N, NADIM
      DOUBLE PRECISION EPSMCH
      DOUBLE PRECISION Z(N), AHESS(NADIM,N), P(N)
C*NS  INTEGER I, IB, IQ, J, JN1, JP1, K
C*NS  INTEGER I, IB,     J, JN1, JP1, K
      DOUBLE PRECISION BETA, DJ, G, GAMMA, GAMMA1, PJ, T
      GAMMA=0.0D+0
      J=1
      DO 30 I=1,N
      T=Z(I)
      IF(I.EQ.1) GOTO 20
      K=I-1
      DO 10 IB=1,K
      T=T-P(IB)*AHESS(I,IB)
 10   CONTINUE
 20   P(I)=T
      GAMMA=GAMMA+T*T/AHESS(I,I)
 30   CONTINUE
      GAMMA1=1.0D+0-GAMMA
      GAMMA=EPSMCH
      IF(GAMMA1.GT.EPSMCH) GAMMA=GAMMA1
      IF(-GAMMA1.GT.EPSMCH) GAMMA=-GAMMA1
      JN1=N+1
      DO 50 I=1,N
      J=JN1-I
      PJ=P(J)
      DJ=AHESS(J,J)
      T=PJ/DJ
      Z(J)=PJ
      BETA=-T/GAMMA
      G=GAMMA+PJ*T
      AHESS(J,J)=DJ*GAMMA/G
      GAMMA=G
      IF(J.EQ.N) GOTO 50
      JP1=J+1
      DO 40 IB=JP1,N
      T=AHESS(IB,J)
      AHESS(IB,J)=T+BETA*Z(IB)
      Z(IB)=Z(IB)+PJ*T
 40   CONTINUE
 50   CONTINUE
      RETURN
      END
