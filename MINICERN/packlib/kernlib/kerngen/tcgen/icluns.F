*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICLUNS (CHV,JL,JR)
C
C CERN PROGLIB# M432    ICLUNS          .VERSION KERNFOR  4.21  890323
C ORIG. 16/01/89, JZ
C
C-    Locate the first 'unseen' character in CHV(JL:JR)

      CHARACTER    CHV*(*)
#include "kerngen/qnatch.inc"
* Ignoring t=pass

      JX = JL
   12 IF (JX.GT.JR)          GO TO 90
#if defined(CERNLIB_QASCII)
      IVAL = ICHAR(CHV(JX:JX))
      IF (IVAL.LT.32)        GO TO 91
      IF (IVAL.GE.127)       GO TO 91
#endif
#if defined(CERNLIB_QEBCDIC)
      IVAL = ICHAR(CHV(JX:JX))
      IF (NATCH(IVAL+1).EQ.0)   GO TO 91
#endif
      JX = JX + 1
      GO TO 12

   90 JX = 0
   91 ICLUNS = JX
      RETURN
      END
