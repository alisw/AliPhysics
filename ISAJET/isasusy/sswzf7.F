#include "isajet/pilot.h"
        REAL FUNCTION SSWZF7(SS)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Drees' function for charged Higgs/sfermion interference
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL SS
      DOUBLE PRECISION S,MS,MW,MZ,MH,AL,BE,SN,RES,
     $EQ,Q,MUS,XL,TERM
C
      S=SS
      MW=TMP(1)
      MZ=TMP(2)
      MH=TMP(3)
      MS=TMP(4)
      AL=TMP(5)
      BE=TMP(6)
      SN=TMP(7)
C
      EQ=(S+MW**2-MZ**2)/2.D0/MW
      Q=SQRT(MAX(0.D0,EQ**2-S))
      MUS=S+MS**2-MZ**2
      XL=LOG((MW*(EQ+Q)-MUS)/(MW*(EQ-Q)-MUS))
      TERM=BE*S*MS**2+SN*AL*MW*MZ*S
      RES=(S*Q*BE/2.D0+TERM*XL/4.D0/MW)/(S-MH**2)
      SSWZF7=RES
      RETURN
      END
