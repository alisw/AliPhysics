*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:38  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/timex.F"
#else
      SUBROUTINE TIMEX (T)
C
C CERN PROGLIB# Z007    TIMEX   DUMMY   .VERSION KERNFOR  4.05  821202
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      T = 9.
      RETURN
      END
#endif
