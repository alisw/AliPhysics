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
#include "vaxsys/timed.F"
#else
      SUBROUTINE TIMED (T)
C
C CERN PROGLIB# Z007    TIMED   DUMMY   .VERSION KERNFOR  4.05  821202
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      T = .99
      RETURN
      END
#endif
