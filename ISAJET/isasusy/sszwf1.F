#include "isajet/pilot.h"
        REAL FUNCTION SSZWF1(E)
C-----------------------------------------------------------------------
C          SSWZBF: ziss -> wiss f fbar
C          Baer's Z2WFUN
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL MW,MZ,M,F1,F2,E,X,Y,WWID
C
      MW=TMP(1)
      MZ=TMP(2)
      X=TMP(3)
      Y=TMP(4)
      M=AMW
      WWID=GAMW
C
      F1=SQRT(MAX(0.,E**2-MW**2))/
     $((MZ**2+MW**2-M**2-2*MZ*E)**2+WWID**2*M**2)
      F2=(X**2+Y**2)*(3*(MZ**2+MW**2)*MZ*E-4*MZ**2*E*E-2*MZ**2*MW**2)
     $-3*(X**2-Y**2)*MZ*MW*(MZ**2+MW**2-2*MZ*E)
      SSZWF1=F1*F2
      RETURN
      END
