#include "isajet/pilot.h"
        REAL FUNCTION SSWZF5(SS)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Baer's XI2FUN
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C
      REAL MW,PI,SS
      DOUBLE PRECISION M1,M2,M3,EQ,Q,XMUS,D,XLOG,S
      DATA PI/3.14159265/,MW/80./
C
      S=SS
      M1=TMP(1)
      M2=TMP(2)
      M3=TMP(3)
      MW=AMW
C
      EQ=(S+M1**2-M3**2)/2./M1
      Q=DSQRT(MAX(0.D0,EQ**2-S))
      XMUS=M2**2+S-M3**2
      D=(M1*(EQ+Q)-XMUS)/(M1*(EQ-Q)-XMUS)
      XLOG=DLOG(D)
      SSWZF5=PI**2/2./M1*M3*S/4./(SS-MW**2)*XLOG
      RETURN
      END
