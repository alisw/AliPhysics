#include "isajet/pilot.h"
        DOUBLE PRECISION FUNCTION SSDLAM(A,B,C)
C-----------------------------------------------------------------------
C          Kinematic function
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
        DOUBLE PRECISION A,B,C
C            Rewrite SSDLAM=A**2+B**2+C**2-2*A*B-2*A*C-2*B*C
        IF(A.GE.B.AND.A.GE.C) THEN
          SSDLAM=(A-B-C)**2-4*B*C
        ELSEIF(B.GE.A.AND.B.GE.C) THEN
          SSDLAM=(B-A-C)**2-4*A*C
        ELSE
          SSDLAM=(C-A-B)**2-4*A*B
        ENDIF
        RETURN
        END
