#include "isajet/pilot.h"
        REAL FUNCTION SSWZF6(EE)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Drees' function for charged Higgs exchange
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL EE
      DOUBLE PRECISION E,MW,MZ,MH,AL,BE,SN,RES
C
      E=EE
      MW=TMP(1)
      MH=TMP(2)
      MZ=TMP(3)
      AL=TMP(4)
      BE=TMP(5)
      SN=TMP(6)
C
      RES=SQRT(E**2-MZ**2)*(MW**2+MZ**2-2*MW*E)*
     ,(E*(AL**2+BE**2)+2*SN*MZ*AL*BE)/
     ,(MW**2+MZ**2-2*MW*E-MH**2)**2
      SSWZF6=RES
      RETURN
      END
