*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMMPW)
#include "mpwgs/datime.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/datime.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/datime.F"
#else
      SUBROUTINE DATIME (ID,IT)
C
C CERN PROGLIB# Z007    DATIME  DUMMY   .VERSION KERNFOR  4.22  890913
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      COMMON /SLATE/ ISLATE(40)

      DO 24 J=1,6
   24 ISLATE(J) = 0

      ID = 790929
      IT = 1200
      RETURN
      END
#endif
