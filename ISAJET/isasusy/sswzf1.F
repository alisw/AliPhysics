#include "isajet/pilot.h"
        REAL FUNCTION SSWZF1(EE)
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
      REAL EE
      DOUBLE PRECISION C1,C2,MWI,MZI,MW,PROP,T1,T2,E
C
      E=EE
      C1=TMP(1)
      C2=TMP(2)
      MWI=TMP(3)
      MZI=TMP(4)
      MW=AMW
C
      PROP=(MWI**2+MZI**2-2*MWI*E-MW**2)**2
      T1=C1*(3*(MWI**2+MZI**2)*MWI*E-4*MWI**2*E*E-2*MWI**2*MZI**2)
      T2=-3*C2*MWI*MZI*(MWI**2+MZI**2-2*MWI*E)
      SSWZF1=SQRT(MAX(0.D0,E*E-MZI**2))/PROP*(T1+T2)
      RETURN
      END
