*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:46  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION ICLOC (CHI,NI,CHV,JL,JR)
C
C CERN PROGLIB# M432    ICLOC           .VERSION KERNFOR  4.21  890323
C ORIG. 09/02/89, JZ
C
C-    Locate CHI of NI characters inside CHV(JL:JR), case sensitive

      CHARACTER    CHI*(*), CHV*(*), CHWK1*1

      JVV   = JL - 1
      JVVE  = JR+1 - NI
      CHWK1 = CHI(1:1)

C--          Scan CHV to find the first char. of CHI

   12 JVV = JVV + 1
      IF (JVV.GT.JVVE)                GO TO 90
      IF (CHV(JVV:JVV).NE.CHWK1)      GO TO 12

C--          Match the further characters

      JVX = JVV
      JI  = 1
   14 JI  = JI + 1
      IF (JI.GT.NI)                   GO TO 91
      JVX = JVX + 1
      IF (CHV(JVX:JVX).EQ.CHI(JI:JI)) GO TO 14
      GO TO 12

   90 JVV = 0
   91 ICLOC = JVV
      RETURN
      END
