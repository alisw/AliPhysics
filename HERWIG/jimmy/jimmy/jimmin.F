

*CMZ :          07/05/95  19.40.21  by  Jonathan Butterworth
*-- Author :    Jonathan Butterworth   08/01/95
      SUBROUTINE JIMMIN
C --------------------------------------------------------------
C Initialise default values for JIMMY. Call just after HWIGIN.
C --------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"

      INTEGER I

 8900 FORMAT(A)

C     Set defaults
      PHAD=300.0

C     Minimum PT of the secondary scatters
      PTJIM=2.5

C     Number of attempts to regenerate a scatter which
C     is initially vetoed by energy conservation 
      MAXMSTRY=100

C     Master flag for MI
      MSFLAG=1

      DO I=1,264
        JMRAD(I)=0.D0
      ENDDO
C     Radii for particles.
C     (zero if not known - i.e. most of them,
C      e+ e- get the photon radius.)
C     Photon
      JMRAD(59)=4.7D-1
C     Proton and Neutron and antis
      JMRAD(73)=7.1D-1
      JMRAD(75)=7.1D-1
      JMRAD(91)=7.1D-1
      JMRAD(93)=7.1D-1
C     e+e-
      JMRAD(121)=JMRAD(59)
      JMRAD(127)=JMRAD(59)

      ANOMOFF=.TRUE.
      JMUEO=1

      RETURN
      END



