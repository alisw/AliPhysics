#include "isajet/pilot.h"
      REAL FUNCTION SSXLAM(A,B,C)
C-----------------------------------------------------------------------
C          Kinematic function
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL A,B,C
C          Rewrite SSXLAM=A**2+B**2+C**2-2*A*B-2*A*C-2*B*C
      IF(A.GE.B.AND.A.GE.C) THEN
        SSXLAM=(A-B-C)**2-4*B*C
      ELSEIF(B.GE.A.AND.B.GE.C) THEN
        SSXLAM=(B-A-C)**2-4*A*C
      ELSE
        SSXLAM=(C-A-B)**2-4*A*B
      ENDIF
      RETURN
      END
