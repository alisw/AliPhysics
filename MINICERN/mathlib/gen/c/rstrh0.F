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
      FUNCTION RSTRH0(X)
#include "gen/def64.inc"
     +      DSTRH0,DSTRH1,D

      SROUND(D)=D+(D-SNGL(D))

      RSTRH0=SROUND(DSTRH0(DBLE(X)))
      RETURN

      ENTRY STRH0(X)
      STRH0=SROUND(DSTRH0(DBLE(X)))
      RETURN

      ENTRY RSTRH1(X)

      RSTRH1=SROUND(DSTRH1(DBLE(X)))
      RETURN

      ENTRY STRH1(X)

      STRH1=SROUND(DSTRH1(DBLE(X)))
      RETURN

      END
#endif
