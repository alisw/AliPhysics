*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DJMNB(AJ,AM,AN,BETA)
      DJMNB = RDJMNB(AJ,AM,AN,BETA)
      END
#endif
