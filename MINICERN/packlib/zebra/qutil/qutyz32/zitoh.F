      SUBROUTINE ZITOH (INTV,HOLL,NP)

#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      INTEGER      INTV(99), HOLL(99), NP(9)

#include "zebra/q_and.inc"
#include "zebra/q_or.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_shiftr.inc"

      N = NP(1)

      DO 39  JWH=1,N
      INTW = INTV(JWH)
      MWH  = 0

      DO 29  JL=1,4
      MWH  = ISHFTL (MWH,8)
      JV   = IAND (INTW,63)
      MWH  = IOR (MWH, NQHOLL(JV))
   29 INTW = ISHFTR (INTW,6)

   39 HOLL(JWH) = MWH
      RETURN
      END
