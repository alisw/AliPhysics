*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUFIND (IT,ICH,JL,JR)
C
C CERN PROGLIB# V304    IUFIND          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION ICH(9)
C
      DO 12 J= JL,JR
      IF (ICH(J).EQ.IT)   GO TO 24
   12 CONTINUE
      J= JR+1
C
   24 IUFIND= J
      RETURN
      END
