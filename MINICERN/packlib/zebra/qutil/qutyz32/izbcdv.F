      SUBROUTINE IZBCDV (NP)

#include "zebra/zceta.inc"
#include "zebra/zkrakc.inc"

      DIMENSION    NP(9)

#include "zebra/q_shiftr.inc"

      N = NP(1)

      DO 29  JWH=1,N
      JV = ISHFTR (IQHOLK(JWH),24)
   29 IQCETK(JWH) = IQCETA(JV+1)
      RETURN
      END
