      SUBROUTINE IZBCDV (NP)

#include "zebra/zceta.inc"
#include "zebra/zkrakc.inc"

#include "zebra/q_and.inc"

      DO 29  JWH=1,NP
      JV = IAND (IQHOLK(JWH),255)
   29 IQCETK(JWH) = IQCETA(JV+1)
      RETURN
      END
