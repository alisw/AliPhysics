*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:39  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:04  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      FUNCTION BETA(X,Y)
      BETA=GAMMA(X)*GAMMA(Y)/GAMMA(X+Y)
      RETURN
      END
