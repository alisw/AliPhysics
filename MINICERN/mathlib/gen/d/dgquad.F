*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
*         Double Prec. version only on 32 bit machines.
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DGQUAD(F,A,B,N)
#include "gen/imp64.inc"
      DIMENSION X(1),W(1)
      EXTERNAL F

      CALL D107D1(1,F,A,B,N,X,W)
      DGQUAD=X(1)
      RETURN
      END
#endif
