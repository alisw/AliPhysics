*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:06  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RSININ(X)
#include "gen/def64.inc"
     +      DSININ,DCOSIN,D

      SROUND(D)=D+(D-SNGL(D))

      RSININ=SROUND(DSININ(DBLE(X)))
      RETURN

      ENTRY SININT(X)
      SININT=SROUND(DSININ(DBLE(X)))
      RETURN

      ENTRY RCOSIN(X)
      RCOSIN=SROUND(DCOSIN(DBLE(X)))
      RETURN

      ENTRY COSINT(X)
      COSINT=SROUND(DCOSIN(DBLE(X)))
      RETURN

      END
#endif
