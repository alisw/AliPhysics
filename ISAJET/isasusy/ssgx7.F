#include "isajet/pilot.h"
        REAL FUNCTION SSGX7(ET)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's XT7 - X- eq. a.6.g of prd45,142 (1992)
C          Modified for t_1 and t_2 eigenstates
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL ET
        DOUBLE PRECISION DET,DMG,DMT,DMZ,DMT1,DMT2,TOP
        DOUBLE PRECISION BOT,PT,DXT7,SSDLAM,PI,TOPS
        DATA PI/3.14159265D0/
        DET=ET
        DMG=TMP(1)
        DMT=TMP(2)
        DMZ=TMP(3)
        DMT1=TMP(4)
        DMT2=TMP(5)
        PT=DSQRT(DET**2-DMT**2)
        TOPS=SSDLAM(DMG**2+DMT**2-2*DMG*DET,DMT**2,DMZ**2)
        TOP=DSQRT(DMAX1(0.D0,TOPS))
        BOT=(DMG**2+DMT**2-2*DMG*DET-DMT1**2)*
     $      (DMG**2+DMT**2-2*DMG*DET-DMT2**2)
        DXT7=PI**2/2.D0*PT*(DMG**2-DMZ**2-2*DMG*DET)/
     $      (DMG**2-2*DMG*DET+DMT**2)*TOP/BOT
        SSGX7=DXT7
        RETURN
        END
