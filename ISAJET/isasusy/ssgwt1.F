#include "isajet/pilot.h"
        REAL FUNCTION SSGWT1(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + tp + bb
C          Baer's FTBW1
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL E
        DOUBLE PRECISION MWI,MG,MT,DFTBW,PT,ET,MST1,MST2,TOP,BOT
C
        ET=E
        MWI=TMP(1)
        MG=TMP(2)
        MT=TMP(3)
        MST1=TMP(6)
        MST2=TMP(7)
C
C            Rewrite PT=DSQRT(ET**2-MT**2)
        PT=DSQRT((ET-MT)*(ET+MT))
        TOP=(MG**2+MT**2-2*MG*ET-MWI**2)**2*ET*PT
        BOT=(MG**2+MT**2-2*MG*ET-MST1**2)*(MG**2+MT**2-2*MG*ET-
     $      MST2**2)*(MG**2+MT**2-2*ET*MG)
        DFTBW=MG*TOP/BOT
        SSGWT1=DFTBW
        RETURN
        END
