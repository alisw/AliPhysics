*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:47  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/incbyt.F"
#else
      FUNCTION INCBYT (INC,A,JTH,NBITS)
C
C CERN PROGLIB# M423    INCBYT          .VERSION KERNFOR  4.07  830624
C ORIG. 10/02/71 CL
C
      DIMENSION A(99),NBITS(3)
      INTEGER   A
C
#include "kerngen/q_jbit.inc"
C
C
      JWORD= (JTH-1)/NBITS(2) + 1
      JPOS = JTH - (JWORD-1)*NBITS(2)
      LOC  = (JPOS-1)*NBITS(1) + 1
      ICT  = INC + JBYT (A(JWORD),LOC,NBITS(1))
      IF (ICT .GT. NBITS(3))     GO TO 2
      CALL SBYT (ICT,A(JWORD),LOC,NBITS)
      INCBYT= 0
      RETURN
C
    2 CALL SBYT (NBITS(3),A(JWORD),LOC,NBITS(1))
      INCBYT= ICT - NBITS(3)
      RETURN
      END
#endif
