*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIRTD)
#include "irtdgs/ie3fod.F"
#elif defined(CERNLIB_QMVAX)&&defined(CERNLIB_QMALPH)
#include "vaxgs/ie3foda.F"
#elif defined(CERNLIB_QMVAX)&&(!defined(CERNLIB_QMALPH))
#include "vaxgs/ie3fodv.F"
#else
      SUBROUTINE IE3FOD (MS,MT,NDPN,JBAD)
C
C CERN PROGLIB# M220    IE3FOD          .VERSION KERNFOR  4.15  861204

C-    Convert double precision for output with copy
C-    from source in native to target in IEEE data format
C-    Dummy routine valid for machines with IEEE internal repr.

      DIMENSION    MS(99), MT(99)

      CALL UCOPY (MS,MT,2*NDPN)

      JBAD = 0
      RETURN
      END
#endif
