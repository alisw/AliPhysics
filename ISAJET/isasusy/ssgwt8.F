#include "isajet/pilot.h"
        REAL FUNCTION SSGWT8(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + tp + bb
C          Baer's FTBW13/  Drees' G_8
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL E
        DOUBLE PRECISION ET,MWI,MG,MT,MST1,MST2,PT,EBMX,EBMN
        DOUBLE PRECISION TOP,BOT,DFTBW
        DOUBLE PRECISION MB,Z,RDL,DEN,R1,R2,R3,SSDLAM
C
        ET=E
        MWI=TMP(1)
        MG=TMP(2)
        MT=TMP(3)
        MST1=TMP(6)
        MST2=TMP(7)
        MB=AMBT
C
C            Rewrite PT=DSQRT(ET**2-MT**2)
        PT=DSQRT((ET-MT)*(ET+MT))
        Z=(MG**2+MT**2-2*MG*ET+MB**2-MWI**2)/2.
        R1=1.D0+MT**2/MG**2-2.D0*ET/MG
        R2=MB**2/MG**2
        R3=MWI**2/MG**2
        RDL=DSQRT(DMAX1(0.D0,SSDLAM(R1,R2,R3)))
        DEN=MG**2+MT**2-2*ET*MG
        EBMX=(2*Z*(1.D0-ET/MG)+PT*MG*RDL)*MG/2.D0/DEN
        EBMN=(2*Z*(1.D0-ET/MG)-PT*MG*RDL)*MG/2.D0/DEN
        TOP=(MG**2+MT**2-2*MG*ET-MWI**2)*(EBMX-EBMN)
        BOT=(MG**2+MT**2-2*MG*ET-MST1**2)*
     $  (MG**2+MT**2-2*MG*ET-MST2**2)
        DFTBW=-MG*MT*TOP/BOT
        SSGWT8=DFTBW
        RETURN
        END
