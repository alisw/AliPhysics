*
* $Id$
*
* $Log$
* Revision 1.1.2.1  2002/07/11 17:14:48  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:34  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if (defined(CERNLIB_DOUBLE))&&(!defined(CERNLIB_IBM))
      FUNCTION GAMMA(X)
      DOUBLE PRECISION DGAMMA,D

      SROUND(D)=D+(D-SNGL(D))
      GAMMA=SROUND(DGAMMA(DBLE(X)))
      RETURN
      END
#endif
