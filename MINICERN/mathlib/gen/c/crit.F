*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION CRIT(SUM,DEL,ACCS)
      COMPLEX SUM,DEL
      CRIT=REAL(DEL)**2+AIMAG(DEL)**2-ACCS*(REAL(SUM)**2+AIMAG(SUM)**2)
      RETURN
      END
