*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CTRANS (CHOLD,CHNEW,CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CTRANS          .VERSION KERNFOR  4.24  900424
C ORIG. 13/03/90, Jamie Shiers
C
C-    Replace all occurrences of CHOLD in CHV(JLP:JRP) by CHNEW

      CHARACTER*(*) CHV
      CHARACTER*1   CHOLD,CHNEW

      DO 24 J=JLP,JRP
   24 IF(CHV(J:J).EQ.CHOLD) CHV(J:J) = CHNEW
      END
