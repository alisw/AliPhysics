#include "isajet/pilot.h"
      SUBROUTINE ORDECR(IA,IB,N)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        return an ordered array (by size of absolute values)
C-        Warning: input array is destroyed
C-
C-   Inputs  : 
C-   IA(N) = input array
C-   Outputs : 
C-   IB(N) = output ordered array
C-
C-   Created   9-MAY-1988   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER IA(*),IB(*),N,I,J,JSEL
C----------------------------------------------------------------------
      DO 2 I=1,N
        JSEL=0
        IB(I)=0
        DO 1 J=1,N
          IF(IABS(IA(J)).GT.IABS(IB(I))) THEN
            IB(I)=IA(J)
            JSEL=J
          ENDIF
    1   CONTINUE
        IF(JSEL.GT.0) IA(JSEL)=0
    2 CONTINUE
  999 RETURN
      END
