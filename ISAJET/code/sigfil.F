#include "isajet/pilot.h"
      SUBROUTINE SIGFIL(SIG,I1,I2,I3,I4)
C          Fill /JETSIG/ arrays if SIG > 0
C          Write error message if SIG < 0
#include "isajet/itapes.inc"
#include "isajet/jetsig.inc"
C
      IF(SIG.GT.0) THEN
        NSIGS=NSIGS+1
        SIGMA=SIGMA+SIG
        SIGS(NSIGS)=SIG
        INOUT(NSIGS)=I1+IOPAK*(I2+IOPAK*(I3+IOPAK*I4))
      ELSEIF(SIG.LT.0.) THEN
        WRITE(ITLIS,1010) SIG,I1,I2,I3,I4
1010    FORMAT(' ERROR IN SIGFIL ... SIG = ',E12.5,' FOR ',4I6)
      ENDIF
      RETURN
      END
