*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:40  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EISST4(NM,N,AR,RLB,RUB,MM,M,WR,IERR,WORK)
C     SOME EIGENVALUES OF A REAL SYMMETRIC TRIDIAGONAL MATRIX
C
      DIMENSION AR(NM,2),WR(MM),WORK(1)
      WORK(1)=AR(1,2)
      IF(N.LT.2) GO TO 150
      DO 100 I=2,N
      WORK(I)=AR(I,2)
      WORK(N+I)=AR(I,1)
  100 WORK(2*N+I)=AR(I,1)**2
  150 EPS1=0.0
      CALL BISECT(N,EPS1,WORK,WORK(N+1),WORK(2*N+1),RLB,RUB,MM,M,WR,
     X WORK(5*N+1),IERR,WORK(3*N+1),WORK(4*N+1))
      RETURN
      END
