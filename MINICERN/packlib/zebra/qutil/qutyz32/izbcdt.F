      SUBROUTINE IZBCDT (NP,ITABT)

#include "zebra/quest.inc"
#include "zebra/zceta.inc"
#include "zebra/zkrakc.inc"
C--------------    END CDE                             --------------
      DIMENSION    NP(9), ITABT(99)

#include "zebra/q_shiftr.inc"

      N     = NP(1)
      LIM   = ITABT(1)
      JGOOD = 0
      JBAD  = 0

      DO 29  JWH=1,N
      JV = ISHFTR (IQHOLK(JWH),24)
      JV = IQTCET(JV+1)
      IF (JV.GT.LIM)               GO TO 27
      JV = ITABT(JV+1)
      IF   (JV+1)            29, 27, 24
   24 JGOOD = JGOOD + 1
      IQCETK(JGOOD) = JV
      GO TO 29

   27 JBAD = JBAD + 1
   29 CONTINUE
      IQUEST(1) = JGOOD
      IQUEST(2) = JBAD
      RETURN
      END
