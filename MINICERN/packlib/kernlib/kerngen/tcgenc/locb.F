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
#include "mpwgs/locb.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/locb.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/locb.F"
#else
      FUNCTION LOCB (IVAR)
C
C CERN PROGLIB# N100    LOCB            .VERSION KERNFOR  4.23  891215
C
C-    This is a default which works on several machines
C
      DIMENSION    IVAR(9)

      LOCB = LOC (IVAR)
      RETURN
      END
#endif
