*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICEQU (CHXA, CHXB, NTX)
C
C CERN PROGLIB# M432    ICEQU           .VERSION KERNFOR  4.39  940228
C ORIG. 08/11/93, JZ
C
C-    Check CHXA and CHXB are identical, return zero if so

      CHARACTER    CHXA(512)*1, CHXB(512)*1


      DO 16  JJ=1,NTX
      IF (CHXA(JJ).NE.CHXB(JJ))    GO TO 19
   16 CONTINUE
      ICEQU = 0
      RETURN

   19 ICEQU = JJ
      RETURN
      END
