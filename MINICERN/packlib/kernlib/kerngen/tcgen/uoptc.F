*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UOPTC (DATA,POSS,IOPTV)
C
C CERN PROGLIB# M502    UOPTC           .VERSION KERNFOR  4.10  850320
C ORIG. 07/11/84, JZ
C
      INTEGER      IOPTV(9)
      CHARACTER    *(*)  DATA, POSS


      NP = LEN(POSS)

      DO 24 J=1,NP
      IOPTV(J) = 0
      IF (INDEX(DATA,POSS(J:J)).NE.0)  IOPTV(J)= 1
   24 CONTINUE
      RETURN
      END
