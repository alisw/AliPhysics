*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:30  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_VAX)||defined(CERNLIB_NORD)||defined(CERNLIB_UNIX)
      FUNCTION IILZ(NW,IA,INC)
C
C CERN PROGLIB# F122    IILZ            .VERSION GEN      2.30  890530
C ORIG. 30/05/89        F. Antonelli IBM, FCA
C
      DIMENSION IA(*)
C
      IRES = 0
      IF(INC.GT.0) THEN
        IBEG = 1
        IEND = (NW-1)*INC+1
      ELSE
        IBEG = (1-NW)*INC+1
        IEND = 1
      END IF
      DO 10 J=IBEG, IEND, INC
        IF(IA(J).NE.0) GO TO 777
        IRES = IRES + 1
  10  CONTINUE
C
 777  IILZ = IRES
C
      END
#endif
