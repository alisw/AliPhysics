      SUBROUTINE IZHNUM (HOLL,INTV,NP)

      INTEGER      INTV(99), HOLL(99)

#include "zebra/q_and.inc"

      DO 39  JWH=1,NP
   39 INTV(JWH) = IAND (HOLL(JWH), 255)
      RETURN
      END
