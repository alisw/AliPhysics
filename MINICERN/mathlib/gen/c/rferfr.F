*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:02  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RFERDR(X,K)
#include "gen/def64.inc"
     +    DFERDR,D

      SROUND(D)=D+(D-SNGL(D))

      RFERDR=SROUND(DFERDR(DBLE(X),K))
      RETURN
 
      ENTRY FERDR(X,K)
      FERDR=SROUND(DFERDR(DBLE(X),K))
      RETURN

      END
#endif
