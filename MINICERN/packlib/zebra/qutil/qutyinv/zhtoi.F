      SUBROUTINE ZHTOI (HOLL,INTV,NP)

#include "zebra/zceta.inc"
C--------------    END CDE                             --------------
      INTEGER      INTV(99), HOLL(99)

#include "zebra/q_and.inc"
#include "zebra/q_or.inc"

      DO 39  JWH=1,NP
      MWH  = HOLL(JWH)
      INTW = 0

      DO 29  JL=1,4
      INTW = ISHFT (INTW,-6)
      JV   = IAND (MWH,255)
      IF (JV.EQ.32)  THEN
          IF (JL.NE.1)             GO TO 29
        ENDIF
      JV   = IQTCET(JV+1)
      INTW = IOR (INTW,ISHFT(JV,18))
   29 MWH  = ISHFT (MWH,-8)

   39 INTV(JWH) = INTW
      RETURN
      END
