*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      REAL FUNCTION FUN(N,RX)
      DOUBLE PRECISION XDB(15), DFUN, FUNX
      REAL RX(N),SNGL
      EXTERNAL DFUN
      DO 10 I=1,N
      XDB(I)=DBLE(RX(I))
 10   CONTINUE
      FUNX=DFUN(N,XDB)
      FUN=SNGL(FUNX)
      RETURN
      END
