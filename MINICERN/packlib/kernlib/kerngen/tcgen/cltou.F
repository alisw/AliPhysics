*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CLTOU (CHV)
C
C CERN PROGLIB# M432    CLTOU           .VERSION KERNFOR  4.21  890323
C ORIG. 11/02/86 A. PETRILLI
C NEW    9/02/89 JZ, for speed
C
C-    Convert character string CHV from lower to upper case.

      CHARACTER    CHV*(*)

#include "kerngen/qnatch.inc"
* Ignoring t=pass


      DO 19  JJ=1,LEN(CHV)
#if defined(CERNLIB_QASCII)
          J = ICHAR(CHV(JJ:JJ))
          IF (J.LT.97)       GO TO 19
          IF (J.GE.123)      GO TO 19
          CHV(JJ:JJ) = CHAR(J-32)
#endif
#if defined(CERNLIB_QEBCDIC)
          J  = ICHAR(CHV(JJ:JJ))
          IF (NATCH(J+1).NE.3)  GO TO 19
          CHV(JJ:JJ) = CHAR(J+64)
#endif
   19 CONTINUE
      END
