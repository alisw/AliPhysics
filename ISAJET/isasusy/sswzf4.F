#include "isajet/pilot.h"
        REAL FUNCTION SSWZF4(SS)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Baer's XI1FUN
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
C-----------------------------------------------------------------------
      REAL SS,PI
      DOUBLE PRECISION M1,M2,M3,EQ,Q,XMUS,XLOG,TERM,S,MW
      DATA PI/3.14159265/
C
      S=SS
      M1=TMP(1)
      M2=TMP(2)
      M3=TMP(3)
      MW=AMW
C
      EQ=(S+M1**2-M3**2)/2./M1
      IF (EQ**2.GE.S) THEN
        Q=DSQRT(EQ**2-S)
      ELSE
        Q=0.D0
      END IF
      XMUS=M2**2+S-M3**2
      XLOG=DLOG((M1*(EQ+Q)-XMUS)/(M1*(EQ-Q)-XMUS))
      TERM=-.5*M1*EQ*Q-.5*(M2**2-M1**2-S)*Q-
     $.25/M1*(M2**2-M3**2)*(M2**2-M1**2)*XLOG
C          SS can stay single precision below
      SSWZF4=PI**2/2./M1/(SS-MW**2)*TERM
      RETURN
      END
