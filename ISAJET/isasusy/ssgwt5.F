#include "isajet/pilot.h"
        REAL FUNCTION SSGWT5(E)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> wiss + tp + bb
C          Baer's FTBW12; EQ. A.3D OF BTW, MODIFIED FOR MB=/0
C          Drees' G_5
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL E
        DOUBLE PRECISION ET,MWI,MG,MT,MSB,MST,PT,EBMX
        DOUBLE PRECISION BOT,TOP,DFTBW,EBMN,XX
        DOUBLE PRECISION MB,Z,RDL,DEN,R1,R2,R3,SSDLAM
C
        ET=E
        MWI=TMP(1)
        MG=TMP(2)
        MT=TMP(3)
        MSB=TMP(4)
        MST=TMP(6)
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
        XX=(MSB**2+2*MG*EBMX-MG**2)/(MSB**2+2*MG*EBMN-MG**2)
        TOP=(MG**2+MT**2-2*MG*ET-MWI**2)*DLOG(XX)
        BOT=MG**2+MT**2-2*MG*ET-MST**2
        DFTBW=-0.5*MT*TOP/BOT
        SSGWT5=DFTBW
        RETURN
        END
