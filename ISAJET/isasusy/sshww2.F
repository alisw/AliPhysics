#include "isajet/pilot.h"
      DOUBLE PRECISION FUNCTION SSHWW2(XX)
C-----------------------------------------------------------------------
C     SSHWW: hi -> z + z*
C     Bisset's FUNZZ
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstmp.inc"
C
      DOUBLE PRECISION XX,EPLN,PROP,TEMP,FN,AAA,MZ,DELTAZ,MH
C
C
      MZ=AMZ
      DELTAZ=GAMZ
      MH=TMP(1)
C
      EPLN=MZ/MH
      PROP=(1.D0-XX)**2 
      PROP=PROP+(EPLN**2)*DELTAZ**2/MH**2
      TEMP=XX**2-12*XX*EPLN**2 
      TEMP=TEMP+8*EPLN**2+12*EPLN**4 
      AAA=XX**2-4*EPLN**2
      IF(AAA.LT.0) THEN
        AAA=0
      ENDIF
      FN=TEMP*SQRT(AAA)/PROP
      SSHWW2=FN
      RETURN
      END
