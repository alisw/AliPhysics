*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LDLSOL (N,NADIM,CHOL,B,X)
      INTEGER N, NADIM
      DOUBLE PRECISION CHOL(NADIM, N)
      DOUBLE PRECISION B(N), X(N)
      INTEGER I, IM1, J, JJ, JP1
      DOUBLE PRECISION SUM
      X(1)=B(1)
      IF(N.EQ.1) GOTO 30
      DO 20 I=2,N
      SUM=B(I)
      IM1=I-1
      DO 10 J=1,IM1
      SUM=SUM-CHOL(I,J)*X(J)
 10   CONTINUE
      X(I)=SUM
 20   CONTINUE
 30   IF(CHOL(N,N).LE.0.0D+0) RETURN
      X(N)=X(N)/CHOL(N,N)
      IF(N.EQ.1) RETURN
      DO 50 JJ=2,N
      J=N-JJ+1
      IF(CHOL(J,J).LE.0.0D+0) RETURN
      SUM=X(J)/CHOL(J,J)
      JP1=J+1
      DO 40 I=JP1,N
      SUM=SUM-CHOL(I,J)*X(I)
 40   CONTINUE
      X(J)=SUM
 50   CONTINUE
      RETURN
      END
