#include "isajet/pilot.h"
      SUBROUTINE SSSAVE(IIN,GAM,IOUT1,IOUT2,IOUT3,IOUT4,IOUT5)
C-----------------------------------------------------------------------
C     Store a SUSY decay mode in /SSMODE/
C     Ver 7.14: Increment NSSMOD only after test
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
C
      INTEGER IIN,IOUT1,IOUT2,IOUT3,IOUT4,IOUT5,I
      REAL GAM
C
      IF (GAM.LE.0.) THEN
        IF(GAM.LT.0.) THEN
          WRITE(LOUT,1000) IIN,IOUT1,IOUT2,IOUT3,IOUT4,IOUT5,GAM
1000      FORMAT(' WARNING: SSSAVE: ',I5,' --> ',5I5,E14.5)
        ENDIF
        LSSMOD=.FALSE.
        GO TO 999
      ENDIF
      NSSMOD=NSSMOD+1
      LSSMOD=.TRUE.
      IF(NSSMOD.GT.MXSS) THEN
        WRITE(LOUT,*) 'SSSAVE: TOO MANY MODES, IIN = ',IIN
        STOP99
      ENDIF
      ISSMOD(NSSMOD)=IIN
      JSSMOD(1,NSSMOD)=IOUT1
      JSSMOD(2,NSSMOD)=IOUT2
      JSSMOD(3,NSSMOD)=IOUT3
      JSSMOD(4,NSSMOD)=IOUT4
      JSSMOD(5,NSSMOD)=IOUT5
      GSSMOD(NSSMOD)=GAM
      BSSMOD(NSSMOD)=0.
      MSSMOD(NSSMOD)=0
C          Check that quarks and gluons appear at end of list.
      DO 100 I=1,4
        IF(IABS(JSSMOD(I,NSSMOD)).LE.9.AND.
     $  IABS(JSSMOD(I+1,NSSMOD)).GT.9) THEN
          WRITE(LOUT,1100) IIN,IOUT1,IOUT2,IOUT3,IOUT4,IOUT5
1100      FORMAT(' WARNING: SSSAVE: BAD ORDER: ',I5,' --> ',5I5)
          STOP99
        ENDIF
100   CONTINUE
C
999   RETURN
      END
