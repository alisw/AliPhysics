      FUNCTION IZBCD (MCHAR)

#include "zebra/zceta.inc"

#include "zebra/q_shiftr.inc"

      JV = ISHFTR (MCHAR,24)
      IZBCD = IQCETA(JV+1)
      RETURN
      END
