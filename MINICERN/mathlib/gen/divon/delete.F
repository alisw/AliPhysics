*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DELETE (NLEN,IARRAY,IVALUE,IPOS)
      INTEGER NLEN, IVALUE, IPOS
      INTEGER IARRAY(NLEN)
      INTEGER I, NLENM1
      IPOS=0
      DO 10 I=1,NLEN
      IF(IVALUE.NE.IARRAY(I)) GOTO 10
      IPOS=I
      GOTO 20
 10   CONTINUE
      RETURN
 20   IF(IPOS.EQ.NLEN) RETURN
      NLENM1=NLEN-1
      DO 30 I=IPOS,NLENM1
      IARRAY(I)=IARRAY(I+1)
 30   CONTINUE
      RETURN
      END
