*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:11  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION CELINT(Z,AKC,A,B)
C
      COMPLEX CELINT,Z
      COMPLEX*16 WELINT,W
      DOUBLE PRECISION D

      SROUND(D)=D+(D-SNGL(D))
      W=WELINT(DCMPLX(DBLE(REAL(Z)),DBLE(AIMAG(Z))),DBLE(AKC),DBLE(A),
     1                DBLE(B))
      CELINT=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))
      RETURN
      END
#endif
