*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE CBSJA(Z,A,NMAX,ND,CB)
      COMPLEX Z,CB(0:*)
#include "gen/defc64.inc"
     + V(0:100),W
#include "gen/def64.inc"
     + D
#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      SROUND(D)=D+(D-SNGL(D))
      W=Z
      CALL WBSJA(W,DBLE(A),NMAX,ND,V)
      DO 1 N = 0,NMAX
    1 CB(N)=CMPLX(SROUND(DREAL(V(N))),SROUND(DIMAG(V(N))))
      RETURN
      END
#endif
