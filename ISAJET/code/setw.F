#include "isajet/pilot.h"
      SUBROUTINE SETW
C
C          Set the W parameters in /WCON/.
C          SIN2W         = sin**2(theta-sub-w)
C          AQ, BQ        = vector, axial couplings normalized to ALFA.
C          MATCH(IQ1,IW) = Cabibbo favored type for W --> QK1 + QK2.
C          WCBR(IQ,IW)   = cumulative branching ratio for JETTYP(1)=IQ
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/wcon.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/nodcay.inc"
#include "isajet/const.inc"
#include "isajet/xmssm.inc"
C
      REAL SINW,COSW,AMW,AMZ,AW,FACZ,GAMW,GAMZ,TERM,SUM,AM1,AMASS,AM2
      INTEGER I1,I2,I3,J,INDEX,IFL,NGAM,NUP,IW,IQ1,IQ2,IFL1,JET,IQ,IFL2
      INTEGER IW1
      REAL T3(12),EQ3(12)
      INTEGER NUTYP(25),LISTJ(25)
#if defined(CERNLIB_SINGLE)
      REAL SIN2WD,SINWD,COSWD,AWD,FACZD
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION SIN2WD,SINWD,COSWD,AWD,FACZD
#endif
      DATA T3/.5,-.5,-.5,.5,-.5,.5,.5,-.5,.5,-.5,.5,-.5/
      DATA EQ3/2.,-1.,-1.,2.,-1.,2.,0.,-3.,0.,-3.,0.,-3./
      DATA NUTYP/13*0,1,1,0,0,1,1,0,0,1,1,0,0/
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
C
C          Masses can be changed with WMASS
C
      SINW=SQRT(SIN2W)
      COSW=SQRT(1.-SIN2W)
      AMW=WMASS(2)
      AMZ=WMASS(4)
C
C          Couplings for Weinberg-Salam model
C
      AW=1./(2.*SQRT2*SINW)
      FACZ=1./(2.*SINW*COSW)
      EZ=SQRT((1.-SIN2W)/SIN2W)
      DO 110 IFL=1,12
        AQ(IFL,1)=EQ3(IFL)/3.
        BQ(IFL,1)=0.
        AQ(IFL,2)=AW
        BQ(IFL,2)=AW
        AQ(IFL,3)=AW
        BQ(IFL,3)=AW
        AQ(IFL,4)=FACZ*(T3(IFL)-2.*EQ3(IFL)/3.*SIN2W)
        BQ(IFL,4)=FACZ*T3(IFL)
110   CONTINUE
#if defined(CERNLIB_SINGLE)
C          Double precision couplings not needed.
      EZDP=EZ
      DO 120 IW=1,4
      DO 120 IFL=1,12
        AQDP(IFL,IW)=AQ(IFL,IW)
        BQDP(IFL,IW)=BQ(IFL,IW)
120   CONTINUE
#endif
#if defined(CERNLIB_DOUBLE)
C          Double precision couplings for 32-bit machines.
      SIN2WD=SIN2W
      SINWD=DSQRT(SIN2WD)
      COSWD=DSQRT(1.-SIN2WD)
      AWD=1./(2.*DSQRT(2.D0)*SINWD)
      FACZD=1./(2.*SINWD*COSWD)
      EZDP=COSWD/SINWD
      DO 120 IFL=1,12
        AQDP(IFL,1)=EQ3(IFL)/3.D0
        BQDP(IFL,1)=0.
        AQDP(IFL,2)=AWD
        BQDP(IFL,2)=AWD
        AQDP(IFL,3)=AWD
        BQDP(IFL,3)=AWD
        AQDP(IFL,4)=FACZD*(T3(IFL)-2.D0*EQ3(IFL)/3.D0*SIN2WD)
        BQDP(IFL,4)=FACZD*T3(IFL)
120   CONTINUE
#endif
C
C          Widths
C
      NGAM=12
      IF(AMLEP(5)+AMLEP(6).GT.AMW) NGAM=9
      GAMW=GF*AMW**3/(6.*PI*SQRT2)*NGAM
      NUP=3
      IF(2.*AMLEP(6).GT.AMZ) NUP=2
      GAMZ=NUP*3.*(AQ(1,4)**2+BQ(1,4)**2)+3.*3.*(AQ(2,4)**2+BQ(2,4)**2)
     1+3.*(AQ(7,4)**2+BQ(7,4)**2+AQ(8,4)**2+BQ(8,4)**2)
      GAMZ=GAMZ*2./FACZ**2
      GAMZ=GAMZ*GF*AMZ**3/(12.*PI*SQRT2)
      WGAM(1)=0.
      WGAM(2)=GAMW
      WGAM(3)=GAMW
      WGAM(4)=GAMZ
C
C          Branching ratios for secondary W+- and Z0
C
      DO 210 IW=2,4
        IW1=IW-1
        SUM=0.
        CUMWBR(1,IW1)=0.
        DO 220 IQ1=2,25
          CUMWBR(IQ1,IW1)=CUMWBR(IQ1-1,IW1)
          IQ2=MATCH(IQ1,IW)
          IF(IQ2.EQ.0) GO TO 220
          IF(.NOT.(GOWMOD(IQ1,IW-1).AND.GOWMOD(IQ2,IW-1))) GO TO 220
          IFL1=LISTJ(IQ1)
          IFL2=LISTJ(IQ2)
          AM1=AMASS(IFL1)
          AM2=AMASS(IFL2)
          IF(AM1+AM2.GE.WMASS(IW)) GO TO 220
          TERM=AQ(IQ1/2,IW)**2+BQ(IQ1/2,IW)**2
          IF(IQ1.LE.13) TERM=3.*TERM
          CUMWBR(IQ1,IW1)=CUMWBR(IQ1-1,IW1)+TERM
          SUM=SUM+TERM
220     CONTINUE
        IF(SUM.LE.0.) THEN
          WRITE(ITLIS,2000) IW
2000      FORMAT(//' ***** NO ALLOWED DECAY MODE FOR SECONDARY W TYPE',
     $    I2,' *****')
          STOP 99
        ENDIF
        DO 230 IQ1=2,25
          CUMWBR(IQ1,IW1)=CUMWBR(IQ1,IW1)/SUM
230     CONTINUE
210   CONTINUE
C
C          Decay channels for DRELLYAN
C
      IF(KEYS(3)) THEN
        DO 310 IW=1,4
          COUT(IW)=0.
          IF(.NOT.GODY(IW)) GO TO 310
          DO 320 IQ1=2,25
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 320
            IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 320
            IF(NUTYP(IQ1)*NUTYP(IQ2).EQ.1.AND.NONUNU) GO TO 320
            IFL1=IQ1/2
            TERM=.5*(AQ(IFL1,IW)**2+BQ(IFL1,IW)**2)
            IF(IQ1.LE.13) TERM=3.*TERM
            COUT(IW)=COUT(IW)+TERM
320       CONTINUE
          IF(COUT(IW).EQ.0.) THEN
            WRITE(ITLIS,3000) IW
3000        FORMAT(//' ***** ERROR IN SETW ... NO ALLOWED DECAY MODE ',
     $      'FOR W TYPE',I2,' *****')
            STOP 99
          ENDIF
310     CONTINUE
C          W branching ratios
        DO 330 IW=1,4
          IF(.NOT.GODY(IW)) GO TO 330
          SUM=0.
          DO 340 IQ1=1,25
            WCBR(IQ1,IW)=SUM
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 340
            IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 340
            IF(NUTYP(IQ1)*NUTYP(IQ2).EQ.1.AND.NONUNU) GO TO 340
            IFL1=IQ1/2
            TERM=.5*(AQ(IFL1,IW)**2+BQ(IFL1,IW)**2)/COUT(IW)
            IF(IQ1.LE.13) TERM=3.*TERM
            SUM=SUM+TERM
            WCBR(IQ1,IW)=SUM
340       CONTINUE
330     CONTINUE
      ENDIF
C
C          Calculate branching ratios for WPAIR events summed over
C          modes allowed by WMODE cards.
C          TBRWW = total allowed branching ratio.
C          RBRWW = relative branching ratios.
C          TBRWW*RBRWW = physical branching ratios.
C
      IF((KEYS(2).AND.(.NOT.GOMSSM)).OR.KEYS(6)
     ,.OR.KEYS(7).OR.KEYS(9).OR.KEYS(10)) THEN
        DO 400 JET=1,2
          TBRWW(1,JET)=1.
          DO 410 IW=2,4
            TBRWW(IW,JET)=0.
            IF(KEYS(6).OR.KEYS(9)) THEN
              IF(.NOT.GOQ(IW,JET)) GO TO 410
            ELSEIF((KEYS(2).OR.KEYS(7).OR.KEYS(10)).AND..NOT.GOMSSM)THEN
              IF(.NOT.GOQ(IW+25,JET)) GO TO 410
            ELSEIF((KEYS(7).OR.KEYS(10)).AND.GOMSSM) THEN
              IF(.NOT.GOQ(IW+76,JET)) GO TO 410
            ENDIF
            SUM=0.
            DO 420 IQ=1,12
              RBRWW(IQ,IW,JET)=0.
              IQ1=2*IQ
              IQ2=MATCH(IQ1,IW)
              IF(IQ2.EQ.0) GO TO 420
              IFL1=IQ1/2
              IF(IQ1.GT.13) IFL1=IFL1+4
              IFL2=IQ2/2
              IF(IQ2.GT.13) IFL2=IFL2+4
              AM1=AMASS(IFL1)
              AM2=AMASS(IFL2)
              IF(AM1+AM2.GE.WMASS(IW)) GO TO 420
              TERM=AQ(IQ1/2,IW)**2+BQ(IQ1/2,IW)**2
              IF(IQ1.LE.13) TERM=3*TERM
              SUM=SUM+TERM
              IF(.NOT.(GOWW(IQ1,JET).AND.GOWW(IQ2,JET))) GO TO 420
              RBRWW(IQ,IW,JET)=TERM
              TBRWW(IW,JET)=TBRWW(IW,JET)+TERM
420         CONTINUE
            TBRWW(IW,JET)=TBRWW(IW,JET)/SUM
            IF(TBRWW(IW,JET).GT.0.) THEN
              DO 430 IQ=1,12
430           RBRWW(IQ,IW,JET)=RBRWW(IQ,IW,JET)/(SUM*TBRWW(IW,JET))
            ELSE
              WRITE(ITLIS,445) IW,JET
445           FORMAT(/' ***** NO ALLOWED MODE FOR W TYPE ',I2,
     $        ' IN JET ',I2,' *****'/)
              STOP 99
            ENDIF
410       CONTINUE
400   CONTINUE
      ENDIF
      RETURN
      END
