*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DAWSON(X)

      DAWSON=RDAWSN(X)
      RETURN

      END
#endif
