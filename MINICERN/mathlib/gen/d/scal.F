*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:21  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SCAL (M,ER)
C-----SUBROUTINE SCAL CALCULATES ER=VARIANCE OF THEORETICAL FUNCTION Y
#include "d510pl.inc"
#include "d510si.inc"
      ER=0.
      N=0
      DO 2 J=1,M
      IF (PL(J)) 2,2,1
 1    N=N+1
      DF(N)=DF(J)
 2    CONTINUE
      IF (N.EQ.0) RETURN
      DO 6 I=1,N
      DO 6 L=1,N
      IF (I-L) 3,3,4
 3    K=L*(L-1)/2+I
      GO TO 5
 4    K=I*(I-1)/2+L
 5    ER=ER+Z(K)*DF(I)*DF(L)
 6    CONTINUE
      RETURN
      END
