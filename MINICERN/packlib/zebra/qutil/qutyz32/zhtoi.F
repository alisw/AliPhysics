      SUBROUTINE ZHTOI (HOLL,INTV,NP)

#include "zebra/zceta.inc"
C--------------    END CDE                             --------------
      INTEGER      INTV(99), HOLL(99), NP(9)
#if defined(CERNLIB_QASCII)
      PARAMETER   (IBLANK=32)
#endif
#if defined(CERNLIB_QEBCDIC)
      PARAMETER   (IBLANK=64)
#endif

#include "zebra/q_and.inc"
#include "zebra/q_or.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_shiftr.inc"

      N = NP(1)

      DO 39  JWH=1,N
      MWH  = HOLL(JWH)
      INTW = 0

      DO 29  JL=1,4
      INTW = ISHFTL (INTW,6)
      JV   = IAND (MWH,255)
      IF (JV.EQ.IBLANK)  THEN
          IF (JL.NE.4)             GO TO 29
        ENDIF
      INTW = IOR (INTW, IQTCET(JV+1))
   29 MWH  = ISHFTR (MWH,8)

   39 INTV(JWH) = INTW
      RETURN
      END
