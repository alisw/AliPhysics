      SUBROUTINE IZCHAV (NP)

#include "zebra/zceta.inc"
#include "zebra/zkrakc.inc"
C--------------    END CDE                             --------------
      DIMENSION    NP(9)

#include "zebra/q_and.inc"
#include "zebra/q_shiftr.inc"

      N   = ISHFTR (NP(1)+3, 2)
      JCE = 0

      DO 29  JWH=1,N
      MWH = IQHOLK(JWH)

      JV  = ISHFTR (MWH,24)
      IQCETK(JCE+1) = IQCETA(JV+1)

      JV  = IAND (ISHFTR(MWH,16), 255)
      IQCETK(JCE+2) = IQCETA(JV+1)

      JV  = IAND (ISHFTR(MWH,8),  255)
      IQCETK(JCE+3) = IQCETA(JV+1)

      JV  = IAND (MWH, 255)
      IQCETK(JCE+4) = IQCETA(JV+1)
   29 JCE = JCE + 4
      RETURN
      END
