*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:29  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (defined(CERNLIB_QF2C))&&(!defined(CERNLIB_LINUX))
      DOUBLE PRECISION FUNCTION GREAL(Z)
      DOUBLE PRECISION Z(2)
      GREAL=Z(1)
      RETURN
      END
#endif
