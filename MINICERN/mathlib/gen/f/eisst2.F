*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:40  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EISST2(NM,N,AR,WR,IERR,WORK)
C     ALL EIGENVALUES OF A REAL SYMMETRIC TRIDIAGONAL MATRIX
C
      DIMENSION AR(NM,2),WR(N),WORK(N)
      DO 100 I=1,N
      WR(I)=AR(I,2)
  100 WORK(I)=AR(I,1)
      CALL IMTQL1(N,WR,WORK,IERR)
      RETURN
      END
