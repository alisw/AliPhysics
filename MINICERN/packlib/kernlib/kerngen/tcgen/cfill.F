*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CFILL (CHI, CHV,JL,JR)
C
C CERN PROGLIB# M432    CFILL           .VERSION KERNFOR  4.21  890323
C ORIG. 06/02/89, JZ
C
C-    Fill CHV(JL:JR) with as many copies of CHI as possible

      CHARACTER    CHV*(*), CHI*(*)

      NI = LEN (CHI)
      JJ = JL

   12 II = 1
   14 IF (JJ.GT.JR)          RETURN
      CHV(JJ:JJ) = CHI(II:II)
      JJ = JJ + 1
      IF (II.EQ.NI)          GO TO 12
      II = II + 1
      GO TO 14
      END
