*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION GAPNC(A,X)

      GAPNC=RGAPNC(A,X)
      RETURN

      ENTRY GAGNC(A,X)
      GAGNC=RGAGNC(A,X)
      RETURN

      END
#endif
