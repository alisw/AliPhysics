*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION ASINH(X)
#include "gen/def64.inc"
     + D,DASINH

      SROUND(D)=D+(D-DBLE(SNGL(D)))
      ASINH=SROUND(DASINH(DBLE(X)))
      RETURN
      END
#endif
