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
#include "mpwgs/locf.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/locf.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/locf.F"
#else
      FUNCTION LOCF (IVAR)
C
C CERN PROGLIB# N100    LOCF            .VERSION KERNFOR  4.34  930114
C
C-    This is a default which works on several machines
C
      DIMENSION    IVAR(9)
*    Number of ADdress Units Per Word for Fortran
C                         and its logarithm base 2
      PARAMETER    (NADUPW=4, LADUPW=2)

      J = LOC(IVAR)
      LOCF = ISHFT (J, -LADUPW)
      END
#endif
