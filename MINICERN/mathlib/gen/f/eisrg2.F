*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:39  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE EISRG2(NM,N,AR,WR,WI,IERR)
C     ALL EIGENVALUES OF A REAL GENERAL MATRIX
C
      DIMENSION AR(NM,NM),WR(N),WI(N)
      CALL BALANC(NM,N,AR,LOW,IGH,WR)
      CALL ELMHES(NM,N,LOW,IGH,AR,WR)
      CALL HQR(NM,N,LOW,IGH,AR,WR,WI,IERR)
      RETURN
      END
