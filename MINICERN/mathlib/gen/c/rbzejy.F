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
      SUBROUTINE RBZEJY(A,N,MODE,REL,X)
#include "gen/def64.inc"
     +      DX,D
      DIMENSION X(*),DX(200)

      SROUND(D)=D+(D-SNGL(D))

      ENTRY BZEJY(A,N,MODE,REL,X)

      CALL DBZEJY(DBLE(A),N,MODE,DBLE(REL),DX)
      DO 1 I = 1,N
    1 X(I)=SROUND(DX(I))
      RETURN
      END
#endif
