#include "isajet/pilot.h"
      DOUBLE PRECISION FUNCTION SSALFS(Q2)  
C-----------------------------------------------------------------------
C     Strong coupling formula from page 201 of Barger and Phillips:
C     (using ALQCD4 for 4 flavor Lambda)
C
C     Bisset's STRCPLH
C-----------------------------------------------------------------------
#include "isajet/sssm.inc"
C
      DOUBLE PRECISION Q2,AS,TH5,TH6,PI
      DATA PI/3.14159265D0/
C
      TH5=4*AMBT**2
      TH6=4*AMTP**2
      IF (Q2.LE.TH5)THEN
        AS=12*PI/(25*LOG(Q2/ALQCD4**2))
      ELSE IF(Q2.GT.TH5.AND.Q2.LE.TH6) THEN
        AS=25*LOG(Q2/ALQCD4**2)-2*LOG(Q2/TH5)
        AS=12*PI/AS
      ELSEIF(Q2.GT.TH6)THEN
        AS=25*LOG(Q2/ALQCD4**2)
        AS=AS-2*(LOG(Q2/TH5)+LOG(Q2/TH6))
        AS=12*PI/AS
      ENDIF
      SSALFS=AS
      RETURN
      END
