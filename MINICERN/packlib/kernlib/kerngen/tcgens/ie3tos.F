*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:10  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIRTD)
#include "irtdgs/ie3tos.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/ie3tos.F"
#else
      SUBROUTINE IE3TOS (MS,MT,NWDO,JBAD)
C
C CERN PROGLIB# M220    IE3TOS          .VERSION KERNFOR  4.15  861204

C-    Convert single precision for input with copy
C-    from source in IEEE to target in native data format
C-    Dummy routine valid for machines with IEEE internal repr.

      DIMENSION    MS(99), MT(99)

      CALL UCOPY (MS,MT,NWDO)

      JBAD = 0
      RETURN
      END
#endif
