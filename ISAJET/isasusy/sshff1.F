#include "isajet/pilot.h"
      DOUBLE PRECISION FUNCTION SSHFF1(MH,MF,NUMH)
C-----------------------------------------------------------------------
C     Calculate QCD radiative correction factor, the square brackets
C     in (4.5) of Drees and Hikasa, Phys. Lett. B240, 455 (1990).
C
C     Bisset's QCDRAD (partial)
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
      DOUBLE PRECISION MH,MF
      DOUBLE PRECISION BETA00,LIXX,LI2PP,LI2MM,TEMP,DELTAP,AS,CF
     $,INPOL,PI,ACAP,DELTAH
      DOUBLE PRECISION DDILOG,SSALFS
      INTEGER NUMH
C
      PI=4*ATAN(1.D0)
      BETA00=SQRT(1-4*MF**2/MH**2)
      LIXX=(1-BETA00)/(1+BETA00)
      LI2PP=DDILOG(LIXX)
      LI2MM=DDILOG(-LIXX)
      TEMP=-3*LOG(1/LIXX)*LOG(2/(1+BETA00))
      TEMP=TEMP-2*LOG(BETA00)*LOG(1/LIXX)
      TEMP=TEMP+4*LI2PP+2*LI2MM
      ACAP=(1+BETA00**2)*TEMP
      ACAP=ACAP-3*BETA00*LOG(4/(1-BETA00**2))
      ACAP=ACAP-4*BETA00*LOG(BETA00)       
      IF (NUMH.EQ.3) THEN
         TEMP=19+2*BETA00**2+3*BETA00**4
         TEMP=TEMP*LOG(1/LIXX)/(16*BETA00)
         DELTAP=TEMP+3*(7-BETA00**2)/8
         DELTAP=DELTAP+ACAP/BETA00
       ELSE
         TEMP=3+34*BETA00**2-13*BETA00**4
         TEMP=TEMP*LOG(1/LIXX)/(16*BETA00**3)
         DELTAH=3*(-1+7*BETA00**2)/(8*BETA00**2)
         DELTAH=DELTAH+TEMP+ACAP/BETA00
       ENDIF
       IF (NUMH.EQ.3) THEN
         INPOL=DELTAP+1.5D0*LOG(MH**2/MF**2)
       ELSE
         INPOL=DELTAH+1.5D0*LOG(MH**2/MF**2)
       ENDIF
       AS=SSALFS(MH**2)
       CF=4.D0/3.D0
       SSHFF1=INPOL*CF*AS/PI+1
      RETURN
      END
