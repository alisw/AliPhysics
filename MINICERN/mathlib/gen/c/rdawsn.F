*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RDAWSN(X)
#include "gen/def64.inc"
     +      DDAWSN,D

      SROUND(D)=D+(D-SNGL(D))

      RDAWSN=SROUND(DDAWSN(DBLE(X)))
      RETURN

      ENTRY DAWSON(X)
      DAWSON=SROUND(DDAWSN(DBLE(X)))
      RETURN
      END
#endif
