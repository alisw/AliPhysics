#include "isajet/pilot.h"
        REAL FUNCTION SSSNWS(EE)
C-----------------------------------------------------------------------
C          SSSNWS: sneutrino->stau_1+f+fb' via W*
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL EE
        DOUBLE PRECISION E,M1,M2,MSN,ML1,WID
        E=EE
        M1=TMP(1)
        M2=TMP(2)
        ML1=TMP(3)
        MSN=TMP(4)
        WID=MSN**2*(E*E-ML1*ML1)**1.5/
     $      (MSN**2+ML1**2-2*MSN*E-M1**2)**2
        SSSNWS=WID
        RETURN
        END
