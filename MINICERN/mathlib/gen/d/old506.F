*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE MINNEW
C
      ENTRY MINTS
      ENTRY MINTL
      ENTRY MINTSD
      ENTRY MINTLD
C
      EXTERNAL FCN
C
      CALL MINUIT(FCN,0)
      RETURN
      END
