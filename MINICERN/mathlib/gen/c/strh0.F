*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION STRH0(X)

      STRH0=RSTRH0(X)
      RETURN

      ENTRY STRH1(X)

      STRH1=RSTRH1(X)
      RETURN

      END
#endif
