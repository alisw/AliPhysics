*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:06  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION SININT(X)

      SININT=RSININ(X)
      RETURN

      ENTRY COSINT(X)
      COSINT=RCOSIN(X)
      RETURN

      END
#endif
