*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZLIMI (LUNPAR,ALIMIP)

C-    SET END-OF-TAPE LIMIT

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcf.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNPAR(9), ALIMIP(9)


      LUN     = LUNPAR(1)
      ALIMIT  = ALIMIP(1)

      CALL FZLOC (LUN,0)
      IF (LQFF.EQ.0)               GO TO 99

      IF   (ALIMIT)          21, 24, 31
   21  Q(KQSP+LQFF+39) = 0.
      IQ(KQSP+LQFF+37) = 0
      IQ(KQSP+LQFF+38) = 0
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.0)
     +WRITE (IQLOG,9022) LUN
 9022 FORMAT ('0FZLIMI.  File at LUN=',I3,' has no output limit.')
#endif
      GO TO 99

   24 ALIMIT = REAL (IQ(KQSP+LQFF+20)) / 1000000.
      ALIMIT = ALIMIT + REAL (IQ(KQSP+LQFF+19)) + Q(KQSP+LQFF+39)
      GO TO 34

   31  Q(KQSP+LQFF+39) = ALIMIT
   34 NUM1 = ALIMIT
      FRAC = ALIMIT - REAL(NUM1)
      NUM2 = 1000000. * FRAC
      IQ(KQSP+LQFF+37) = NUM1
      IQ(KQSP+LQFF+38) = NUM2

#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.0)
     +WRITE (IQLOG,9036) LUN,NUM1,NUM2
 9036 FORMAT ('0FZLIMI.  Output limit for file at LUN=',I3,
     F' IS',I7,' Mwords +',I6,' Words.')
#endif
   99 RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
