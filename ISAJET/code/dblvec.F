#include "isajet/pilot.h"
      SUBROUTINE DBLVEC(P,DP)
C
C          Calculate double precision vector DP for 5-vector P.
C          Exact components are 1,2,5 and larger of +,-
C          Ver 6.44: Always use this, even if IF=SINGLE.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL P(5)
      DOUBLE PRECISION DP(5),DPPL,DPMN
      INTEGER K
C
      DO 100 K=1,5
100   DP(K)=P(K)
      IF(DP(4)+ABS(DP(3)).EQ.0.) RETURN
      IF(DP(3).GT.0.) THEN
        DPPL=DP(4)+DP(3)
        DPMN=(DP(1)**2+DP(2)**2+DP(5)**2)/DPPL
      ELSE
        DPMN=DP(4)-DP(3)
        DPPL=(DP(1)**2+DP(2)**2+DP(5)**2)/DPMN
      ENDIF
      DP(3)=0.5D0*(DPPL-DPMN)
      DP(4)=0.5D0*(DPPL+DPMN)
      RETURN
      END
