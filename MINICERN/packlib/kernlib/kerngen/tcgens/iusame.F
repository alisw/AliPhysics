*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUSAME (IA,JST,N,MIN,JSAME)
C
C CERN PROGLIB# M501    IUSAME          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71 CL
C
      DIMENSION IA(*)
C
C
      J  = JST
   21 NSM= 0
      JSM= J
C
   24 NSM= NSM + 1
      IF (J.EQ.N)  GO TO 31
      J = J + 1
      IF (IA(J).EQ.IA(J-1))  GO TO 24
C
      IF (NSM.LT.MIN)  GO TO 21
      GO TO 32
C
   31 IF (NSM.GE.MIN)  GO TO 32
      NSM= 0
      JSM= N + 1
   32 IUSAME=NSM
      JSAME= JSM
      RETURN
      END
