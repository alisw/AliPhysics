*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MODCHL (N,NADIM,AHESS,ALPHA,Z,IFAIL)
      INTEGER N, NADIM, IFAIL
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION AHESS(NADIM,N), Z(N)
      INTEGER I, IB, IP1, J, K
      DOUBLE PRECISION A, BETA, DB, DI, GAMMA, P1, RMAX, T, W
#if defined(CERNLIB_IBM)||defined(CERNLIB_SINGLE)
      RMAX=7.2D+75
#endif
#if (!defined(CERNLIB_IBM))&&(defined(CERNLIB_DOUBLE))
      RMAX= 1.0D+37
#endif
      IFAIL=0
      A=ALPHA
      K=0
      DO 50 I=1,N
      P1=Z(I)
      DI=AHESS(I,I)
      T=A*P1
      DB=DI+T*P1
      AHESS(I,I)=DB
      IF(DB.GE.1.0D+0) GOTO 10
      IF(DB.GT.0.0D+0.AND.DI.LE.RMAX*DB) GOTO 10
      IFAIL=1
      RETURN
 10   GAMMA=DI/DB
      BETA=T/DB
      A=A*GAMMA
      K=K+I
      J=K
      IF(I.EQ.N) GOTO 50
      IP1=I+1
      IF(GAMMA.GE.2.5D-1) GOTO 30
      DO 20 IB=IP1,N
      T=AHESS(IB,I)
      AHESS(IB,I)=T*GAMMA+BETA*Z(IB)
      Z(IB)=Z(IB)-P1*T
 20   CONTINUE
      GOTO 50
 30   DO 40 IB=IP1,N
      T=AHESS(IB,I)
      W=Z(IB)-P1*T
      Z(IB)=W
      AHESS(IB,I)=BETA*W+T
 40   CONTINUE
 50   CONTINUE
      RETURN
      END
