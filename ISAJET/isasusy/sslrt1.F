#include "isajet/pilot.h"
        REAL FUNCTION SSLRT1(SS)
C-----------------------------------------------------------------------
C          SSLRT1: l_R -> l+tau+stau_1
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL SS
        DOUBLE PRECISION ETMX,ETMN,S,MT1,MT,MLR,BEZI,BEZJ,TM,
     ,AI,AJ,BI,BJ,MZI,MZJ,SNZI,SNZJ,XL,BK1,BK2,BK3,BK,WID,SSDLAM
        S=SS
        MT1=AML1SS
        MT=AMTAU
        MLR=TMP(1)
        BEZI=TMP(2)
        BEZJ=TMP(3)
        AI=TMP(4)
        AJ=TMP(5)
        BI=TMP(6)
        BJ=TMP(7)
        MZI=ABS(TMP(8))
        MZJ=ABS(TMP(9))
        SNZI=SIGN(1.0,TMP(8))
        SNZJ=SIGN(1.0,TMP(9))
        TM=SSDLAM(S,MT**2,MT1**2)
        XL=DSQRT(MAX(0.D0,TM))
        ETMN=(S+MT**2-MT1**2-XL*(MLR**2-S)/(MLR**2+S))*(MLR**2+S)/
     ,       (2*S)/(2*MLR)
        ETMX=(S+MT**2-MT1**2+XL*(MLR**2-S)/(MLR**2+S))*(MLR**2+S)/
     ,       (2*S)/(2*MLR)
        BK1=-(ETMX-ETMN)*((ETMX+ETMN)*MLR*S-
     ,       (S+MT**2-MT1**2)*MLR**2)/2.D0
        BK2=(ETMX-ETMN)*((ETMX+ETMN)*MLR-S-MT**2+MT1**2)/2.D0
        BK3=SNZJ*BI*AJ*MZJ+SNZI*BJ*AI*MZI
        BK=BI*BJ*BK1+AI*AJ*MZI*MZJ*SNZI*SNZJ*BK2+BK3*MT*(MLR**2-S)*
     ,     (ETMX-ETMN)/2.D0
        WID=BEZI*BEZJ*BK/(S-MZI**2)/(S-MZJ**2)
        SSLRT1=WID
        RETURN
        END
