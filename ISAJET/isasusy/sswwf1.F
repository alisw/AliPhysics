#include "isajet/pilot.h"
        REAL FUNCTION SSWWF1(EE)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Baer's WIWFUN
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      DOUBLE PRECISION E,MW2,MW1,MZ,SN,XX,YY,T1,T2,T3,T4
      REAL EE
C
      E=EE
      MW2=TMP(1)
      MW1=TMP(2)
      MZ=TMP(3)
      SN=TMP(4)
      XX=TMP(5)
      YY=TMP(6)
C
      T1=DSQRT(E**2-MW1**2)/(MW2**2+MW1**2-MZ**2-2*MW2*E)**2
      T2=3*E*MW2*(MW2**2+MW1**2)-2*MW2**2*MW1**2-4*MW2**2*E*E
      T3=2*E*MW2-MW2**2-MW1**2
      T4=T1*((XX**2+YY**2)*T2-3*SN*(XX**2-YY**2)*MW2*MW1*T3)
      SSWWF1=T4
      RETURN
      END
