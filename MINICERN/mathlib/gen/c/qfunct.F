*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE QFUNCT(XK,XKPRIM,Q)
      XKPRIM=SQRT(1.-XK**2)
      A=MAX(ABS(XK),XKPRIM)
      EPSI=(1.-A)/(2.*A+4.*SQRT(A)+2.)
      E=EPSI**4
      Q=EPSI+EPSI*E*(2.+E*(15.+E*150.))
      IF(XKPRIM.LT.A) Q=EXP(9.869604401089359/LOG(Q))
      RETURN
      END
