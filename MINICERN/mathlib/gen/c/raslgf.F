*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:04  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RASLGF(MODE,X,M,NL,P)
      DIMENSION P(0:*)
#include "gen/def64.inc"
     +     D,DP(0:100)

      SROUND(D)=D+(D-SNGL(D))

      ENTRY ASLGF(MODE,X,M,NL,P)

      CALL DASLGF(MODE,DBLE(X),M,NL,DP)
      DO 1 N = 0,NL
    1 P(N)=SROUND(DP(N))
      RETURN
      END
#endif
