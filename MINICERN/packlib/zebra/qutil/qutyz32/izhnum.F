      SUBROUTINE IZHNUM (HOLL,INTV,NP)

      INTEGER      INTV(99), HOLL(99), NP(9)

#include "zebra/q_shiftr.inc"

      N = NP(1)

      DO 39  JWH=1,N
   39 INTV(JWH) = ISHFTR (HOLL(JWH), 24)
      RETURN
      END
