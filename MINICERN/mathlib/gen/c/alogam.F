*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_IBM))
      FUNCTION ALOGAM(X)
      ALOGAM = ALGAMA(X)
      END
#endif
