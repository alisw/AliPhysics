*-- Author : JMB
      SUBROUTINE JMEFIN
C -----------------------------------------------------------------
C     Adjust HERWIG cross section in light of Multiple Scattering.
C     Must be called after (from?) HWEFIN.
C -----------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"

      WRITE (JMOUT,1)
    1 FORMAT(/10X,'MODIFIED OUTPUT ON ELEMENTARY PROCESS'/)

      IF (JMUEO.EQ.0) THEN

         IF (NLOST.EQ.0) THEN
            WRITE (6,10)
         ELSE
            AVWGT=AVWGT*(TOTSCAT-NLOST)/TOTSCAT
            WRITE (6,11) INT(TOTSCAT),INT(NLOST),1000.*AVWGT

         ENDIF
      ELSE
         WRITE (6,13)
         WRITE (6,12) INT(TOTSCAT)         
      ENDIF
         
 10   FORMAT(10X,
     &     'NO MULTIPLE SCATTERS:NO MODIFICATION NECESSARY')
 11      FORMAT(1P,
     &        10X,'NUMBER OF SCATTERS = ',I20/
     &        10X,'NO.OF MULT.SCATTERS= ',I20/
     &        10X,'FINAL C-S (PB)     =',G12.4)
 12      FORMAT(1P,
     &        10X,'NUMBER OF SCATTERS = ',I20)
 13      FORMAT(1P,
     &        10X,'MULTIPLE SCATTERS USED FOR UNDERLYING EVENT'/
     &        10X,'NO CHANGE TO TOTAL CROSS SECTION.')
      END
