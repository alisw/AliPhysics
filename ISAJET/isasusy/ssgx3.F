#include "isajet/pilot.h"
        REAL FUNCTION SSGX3(ET)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's XT3 - CHI- eq. a.6.c of prd45,142 (1992)
C          Modified for t_1 and t_2 eigenstates
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL ET
        DOUBLE PRECISION DET,DMG,DMT,DMZ,DMT1,DMT2,TOP,BOT
        DOUBLE PRECISION PT,DXT3,SSDLAM,PI,TOPS
        DATA PI/3.14159265D0/
        DET=ET
        DMG=TMP(1)
        DMT=TMP(2)
        DMZ=TMP(3)
        DMT1=TMP(4)
        DMT2=TMP(5)
        TOPS=SSDLAM(DMG**2+DMT**2-2*DMG*DET,DMT**2,DMZ**2)
        TOP=DSQRT(DMAX1(1.D0,TOPS))
        BOT=(DMG**2+DMT**2-2*DMG*DET-DMT1**2)*
     $      (DMG**2+DMT**2-2*DMG*DET-DMT2**2)
        PT=DSQRT(DET**2-DMT**2)
        DXT3=PI**2*DMG*DET*PT*TOP/BOT/(DMG**2-2*DMG*DET+DMT**2)
        SSGX3=DXT3
        RETURN
        END
