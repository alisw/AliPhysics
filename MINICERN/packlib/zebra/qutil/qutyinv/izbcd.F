      FUNCTION IZBCD (MCHAR)

#include "zebra/zceta.inc"

#include "zebra/q_and.inc"

      JV = IAND (MCHAR,255)
      IZBCD = IQCETA(JV+1)
      RETURN
      END
