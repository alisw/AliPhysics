*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:39  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EISST1(NM,N,AR,WR,ZR,IERR,WORK)
C     ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS OF A REAL
C     SYMMETRIC TRIDIAGONAL MATRIX
C
      DIMENSION AR(NM,2),WR(N),ZR(NM,NM),WORK(N)
      DO 100 I=1,N
      DO 50 J=1,N
   50 ZR(I,J)=0.0
      ZR(I,I)=1.0
      WR(I)=AR(I,2)
  100 WORK(I)=AR(I,1)
      CALL IMTQL2(NM,N,WR,WORK,ZR,IERR)
      RETURN
      END
