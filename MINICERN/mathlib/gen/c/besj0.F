*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION BESJ0(X)
#include "gen/def64.inc"
     +     D,DBESJ0,DBESJ1,DBESY0,DBESY1

      SROUND(D)=D+(D-SNGL(D))
      BESJ0=SROUND(DBESJ0(DBLE(X)))
      RETURN

      ENTRY BESJ1(X)
      BESJ1=SROUND(DBESJ1(DBLE(X)))
      RETURN

      ENTRY BESY0(X)
      BESY0=SROUND(DBESY0(DBLE(X)))
      RETURN

      ENTRY BESY1(X)
      BESY1=SROUND(DBESY1(DBLE(X)))
      RETURN
      END
#endif
