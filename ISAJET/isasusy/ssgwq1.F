#include "isajet/pilot.h"
        REAL FUNCTION SSGWQ1(Q)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + qk + qb
C          The function psi of PRD36, 96 (1987); Eq. 3.2
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL Q,PI
        DOUBLE PRECISION MZ,MG,MS,QS,MGS,MSS,MZS,FUN
        PI=4*ATAN(1.)
        MZ=TMP(1)
        MS=TMP(2)
        MG=AMGLSS
        QS=Q*Q
        MGS=MG*MG
        MZS=MZ*MZ
        MSS=MS*MS
        FUN=QS*(MGS-2*MG*Q-MZS)**2/(MGS-2*MG*Q-MSS)**2/(MGS-2*MG*Q)
        SSGWQ1=PI**2*MG*FUN
        RETURN
        END
