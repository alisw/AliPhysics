#include "isajet/pilot.h"
        REAL FUNCTION SSGWT3(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + tp + bb
C          Drees' G_3
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL E
        DOUBLE PRECISION MWI,MG,MT,EB,BOT,MSB
        DOUBLE PRECISION SSDLAM,TERM,SN,DFTBW
C
        EB=E
        MWI=TMP(1)
        MG=TMP(2)
        MT=TMP(3)
        MSB=TMP(4)
        SN=TMP(8)
C
        
        BOT=(MG**2+AMBT**2-2*MG*EB-MSB**2)**2*
     $       (MG**2+AMBT**2-2*EB*MG)
        TERM=SSDLAM((MG**2+AMBT**2-2*MG*EB),MWI**2,MT**2)
        IF(TERM.GT.0.D0) THEN
          DFTBW=4*MG*MWI*MT*SN*EB*EB/BOT*DSQRT(TERM)
        ELSE
          DFTBW=0.D0
        ENDIF
        SSGWT3=DFTBW
        RETURN
        END
