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
#include "vaxsys/intrac.F"
#elif defined(CERNLIB_QMMPW)
#include "mpwgs/intrac.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/intrac.F"
#else
      LOGICAL FUNCTION INTRAC(DUMMY)
C
C CERN PROGLIB# Z044    INTRAC          .VERSION KERNFOR  4.10  850320
C ORIG. 20/02/83  HRR, DUMMY VERSION
C
      INTRAC = .FALSE.
      RETURN
      END
#endif
