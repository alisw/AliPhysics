#include "isajet/pilot.h"
        REAL FUNCTION SSWZF2(QQ)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Baer's PSIFUN
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL QQ
      DOUBLE PRECISION M1,M2,M3,PI,Q
      DATA PI/3.14159265D0/
C
      Q=QQ
      M1=TMP(1)
      M2=TMP(2)
      M3=TMP(3)
C
      SSWZF2=PI**2*M1*Q*Q*(M1**2-2*M1*Q-M3**2)**2/
     $(M1**2-2*M1*Q-M2**2)**2/(M1**2-2*M1*Q)
      RETURN
      END
