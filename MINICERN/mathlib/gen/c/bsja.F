*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE BSJA(X,A,NMAX,ND,B)
      DOUBLE PRECISION D,DB(0:101)
      DIMENSION B(0:*)

      SROUND(D)=D+(D-SNGL(D))
      CALL DBSJA(DBLE(X),DBLE(A),NMAX,ND,DB)
      GO TO 2

      ENTRY BSIA(X,A,NMAX,ND,B)

      CALL DBSIA(DBLE(X),DBLE(A),NMAX,ND,DB)
    2 DO 1 N = 0,ABS(NMAX)
    1 B(N)=SROUND(DB(N))
      RETURN
      END
#endif
