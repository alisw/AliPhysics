*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:50  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UOPT (DATA,POSS,IOPT,NP)
C
C CERN PROGLIB# M502    UOPT            .VERSION KERNFOR  4.18  880425
C ORIG. 01/08/71
C
      COMMON /SLATE/IP(30),ID(10)
      DIMENSION DATA(9), POSS(9), IOPT(9)
#include "kerngen/wordsize.inc"
C
C
      CALL UBLOW (POSS,IP,NP)
      CALL VZERO (IOPT,NP)
      L   = 1
C
   21 CALL UBLOW (DATA(L),ID,NCHAPW)
C
      DO 24 I=1,NCHAPW
      J = IUCOMP (ID(I),IP,NP)
      IF (J.EQ.0)  RETURN
   24 IOPT(J)= 1
      L = L + 1
      GO TO 21
      END
