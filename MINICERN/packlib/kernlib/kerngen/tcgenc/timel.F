*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/timel.F"
#else
      SUBROUTINE TIMEL (T)
C
C CERN PROGLIB# Z007    TIMEL   DUMMY   .VERSION KERNFOR  4.11  850529
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      T = 9999.
      RETURN
      END
#endif
