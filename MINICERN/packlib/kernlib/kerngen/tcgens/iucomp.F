*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUCOMP (ITEXT,IVECT,N)
C
C CERN PROGLIB# V304    IUCOMP          .VERSION KERNFOR  1.0   710701
C ORIG. 01/01/66 JZ
C
      DIMENSION IVECT(9)
C
      IF (N.EQ.0) GO TO 18
      DO 12 J=1,N
      IF (ITEXT.EQ.IVECT(J))  GO TO 24
   12 CONTINUE
   18 J=0
C
   24 IUCOMP= J
      RETURN
      END
