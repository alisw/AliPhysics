*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/locbyt.F"
#else
      FUNCTION LOCBYT (IT,IA,N,NW,LA,NBITS)
C
C CERN PROGLIB# M428    LOCBYT          .VERSION KERNFOR  4.07  830624
C ORIG. 03/05/72
C
      DIMENSION IA(9)
C
#include "kerngen/q_jbit.inc"
* Ignoring t=pass
C
      IF (N.EQ.0)  GO TO 19
      DO 18 J=1,N,NW
      IF (JBYT(IA(J),LA,NBITS).EQ.IT)  GO TO 21
   18 CONTINUE
   19 LOCBYT= 0
      RETURN
C
   21 LOCBYT= J
      RETURN
      END
#endif
