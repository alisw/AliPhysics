#include "isajet/pilot.h"
        REAL FUNCTION SSGWT2(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + tp + bb
C          Drees' G_2
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL E
        DOUBLE PRECISION MWI,MG,MT,DFTBW,EB,TOP,BOT,MSB
        DOUBLE PRECISION SSDLAM,TERM
C
        EB=E
        MWI=TMP(1)
        MG=TMP(2)
        MT=TMP(3)
        MSB=TMP(4)
C
        TOP=EB**2*(MG**2+AMBT**2-2*MG*EB-MWI**2-MT**2)
        BOT=(MG**2+AMBT**2-2*MG*EB-MSB**2)**2*
     $       (MG**2+AMBT**2-2*EB*MG)
        TERM=SSDLAM((MG**2+AMBT**2-2*MG*EB),MWI**2,MT**2)
        IF(TERM.GT.0.D0) THEN
          DFTBW=MG*TOP/BOT*DSQRT(TERM)
        ELSE
          DFTBW=0.D0
        ENDIF
        SSGWT2=DFTBW
        RETURN
        END
