*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IULOOK(NCHAR,IC,JL,JR)
C
C CERN PROGLIB# M410    IULOOK          .VERSION KERNFOR  3.13  810128
C ORIG. 27/05/67
C
      COMMON /SLATE/ MMM,IBUF(39)
C
      DIMENSION IC(80)
      DATA  IBLANK / 4H     /
C
C
      IBUF(1)=IBLANK
      I=0
      DO 24 J=JL,JR
      IF (IC(J).EQ.IBLANK) GO TO 24
      IF(I.GE.NCHAR) GO TO 29
      I=I+1
      IBUF(I)= IC(J)
   24 CONTINUE
C
   29 IF (I.EQ.0) I=1
      CALL UBUNCH (IBUF,MMM,I)
      IULOOK = MMM
      RETURN
      END
