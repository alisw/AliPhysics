*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUCOLA (ITEXT,IVECT,N)
C
C CERN PROGLIB# V304    IUCOLA          .VERSION KERNFOR  4.18  880425
C ORIG. 25/04/88 JZ
C
      DIMENSION IVECT(9)
C
      IF (N.EQ.0)            GO TO 18
      DO 12   JF=N,1,-1
      IF (ITEXT.EQ.IVECT(JF))    GO TO 24
   12 CONTINUE
   18 JF = 0
C
   24 IUCOLA = JF
      RETURN
      END
