      SUBROUTINE ZITOH (INTV,HOLL,NP)

#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      INTEGER      INTV(99), HOLL(99)

#include "zebra/q_and.inc"
#include "zebra/q_or.inc"

      DO 39  JWH=1,NP
      INTW = INTV(JWH)
      MWH  = 0

      DO 29  JL=1,4
      MWH  = ISHFT (MWH,-8)
      JV   = IAND (INTW,63)
      MWH  = IOR (MWH, ISHFT(NQHOLL(JV),24))
   29 INTW = ISHFT (INTW,-6)

   39 HOLL(JWH) = MWH
      RETURN
      END
