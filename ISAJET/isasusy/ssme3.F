#include "isajet/pilot.h"
      SUBROUTINE SSME3(KTYP,AM,ZI,ZF)
C
C          Give matrix element data for mode most recently saved by 
C          SSSAVE. Call this once for each pole in the matrix element,
C          giving the pole type, mass, and couplings. See /DKYSS3/
C          for more comments.
C
C          Assumes SUSY decay product is always FIRST.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/dkyss3.inc"
C
      INTEGER KTYP,I
      REAL AM
      COMPLEX ZI(2),ZF(2)
C
C          If last SSSAVE failed, then skip the matrix element
C
      IF(.NOT.LSSMOD) RETURN
C
C          If MSSMOD(NSSMOD)=0, have not booked any poles yet for
C          last mode saved. Increment mode counter, and set initial and
C          final poles to next one.
C
      IF(MSSMOD(NSSMOD).EQ.0) THEN
        NMSS3=NMSS3+1
        IF(NMSS3.GT.MXMSS3) THEN
          WRITE(LOUT,*) 'ERROR IN SSME3...TOO MANY MODES=',NMSS3
          STOP99
        ENDIF
        MSSMOD(NSSMOD)=-NMSS3
        J1SS3(NMSS3)=NPSS3+1
        J2SS3(NMSS3)=NPSS3+1
        WTSS3(NMSS3)=0
      ENDIF
C
C          Add pole to list and set second counter to last pole
C
      NPSS3=NPSS3+1
      IF(NPSS3.GT.MXPSS3) THEN
        WRITE(LOUT,*) 'ERROR IN SSME3...TOO MANY POLES=',NPSS3
        STOP99
      ENDIF
      KSS3(NPSS3)=KTYP
      AMSS3(NPSS3)=AM
      DO 100 I=1,2
        ZISS3(I,NPSS3)=ZI(I)
        ZFSS3(I,NPSS3)=ZF(I)
100   CONTINUE
      J2SS3(NMSS3)=NPSS3
C
      RETURN
      END
