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
#include "irtdgs/ie3fos.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/ie3fos.F"
#else
      SUBROUTINE IE3FOS (MS,MT,NWDO,JBAD)
C
C CERN PROGLIB# M220    IE3FOS          .VERSION KERNFOR  4.15  861204

C-    Convert single precision for output with copy
C-    from source in native to target in IEEE data format
C-    Dummy routine valid for machines with IEEE internal repr.

      DIMENSION    MS(99), MT(99)

      CALL UCOPY (MS,MT,NWDO)

      JBAD = 0
      RETURN
      END
#endif
