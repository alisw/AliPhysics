#include "isajet/pilot.h"
      FUNCTION ALQCD(Q2)  
C-----------------------------------------------------------------------
C     Strong coupling formula from page 201 of Barger and Phillips:
C     (using ALQCD4 for 4 flavor Lambda)
C-----------------------------------------------------------------------
      REAL Q2,AS,TH5,TH6,PI,ALQCD4
      LOGICAL FIRST
      SAVE FIRST,PI,TH5,TH6,ALQCD4
      DATA FIRST/.TRUE./
C
      IF(FIRST) THEN
        PI=4.*ATAN(1.)
        TH5=4*AMASS(5)**2
        TH6=4*AMASS(6)**2
        ALQCD4=0.177
        FIRST=.FALSE.
      ENDIF
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
      ALQCD=AS
      RETURN
      END
