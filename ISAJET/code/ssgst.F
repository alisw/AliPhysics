#include "isajet/pilot.h"
        REAL FUNCTION SSGST(S,AMSQ,Z,I,J)
C-----------------------------------------------------------------------
C          Function for Sig(qqbar->z_i + z_j
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
        REAL S,AMSQ,K,Z,MZI,MZJ,RS,TP,BT
        INTEGER I,J,ITHI,ITHJ
C
        MZI=ABS(AMZISS(I))
        MZJ=ABS(AMZISS(J))
        IF (AMZISS(I).LT.0.) THEN
          ITHI=1
        ELSE
          ITHI=0
        END IF
        IF (AMZISS(J).LT.0.) THEN
          ITHJ=1
        ELSE
          ITHJ=0
        END IF
        RS=SQRT(S)
        K=SQRT(S*S+(MZI**2-MZJ**2)**2-2*S*(MZI**2+MZJ**2))/
     $   2./RS
        TP=S*S-(MZI**2-MZJ**2)**2-4*K*S**1.5*Z+4*K*K*S*Z*Z+
     $   4*(-1.)**(ITHI+ITHJ+1)*MZI*MZJ*S
        BT=(S-MZI**2-MZJ**2)/2.-RS*K*Z+AMSQ**2
        SSGST=TP/BT
        RETURN
        END
