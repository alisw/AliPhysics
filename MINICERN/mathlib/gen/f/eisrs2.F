*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:39  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EISRS2(NM,N,AR,WR,IERR,WORK)
C     ALL EIGENVALUES OF A REAL SYMMETRIC MATRIX
C
      DIMENSION AR(NM,NM),WR(N),WORK(1)
      CALL TRED1(NM,N,AR,WR,WORK,WORK)
      CALL TQL1(N,WR,WORK,IERR)
      RETURN
      END
