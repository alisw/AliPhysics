*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:02  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION FRSIN(X)
      FRSIN=RFRSIN(X)
      RETURN
      ENTRY FRCOS(X)
      FRCOS=RFRCOS(X)
      END
#endif
