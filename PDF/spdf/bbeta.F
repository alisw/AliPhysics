*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:38  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:04  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
#if defined(CERNLIB_SINGLE)
      FUNCTION BBETA(X,Y)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION FUNCTION BBETA(X,Y)
#endif
#include "pdf/expdp.inc"
     +       DGAMMA,
     +       X,Y
#if defined(CERNLIB_SINGLE)
      BBETA=GAMMA(X)*GAMMA(Y)/GAMMA(X+Y)
#endif
#if defined(CERNLIB_DOUBLE)
      BBETA=DGAMMA(X)*DGAMMA(Y)/DGAMMA(X+Y)
#endif
      RETURN
      END
