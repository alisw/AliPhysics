#include "isajet/pilot.h"
      SUBROUTINE SIGSSE
C
C          Compute d(sigma)/d(cos theta) for
C          e+ e- ----> SUSY particles
C          See Baer et. al., IJMP A4, 4111 (1989) for sigma's
C          Polarized cross sections added 9/18/95 hb
C          Mixed sbottoms and staus included 10/23/96 hb
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE and WTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
C
C          Extra factor of 1/2 needed because all jets are treated
C          as identical.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetsig.inc"
#include "isajet/eepar.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/q1q2.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstype.inc"
#include "isajet/brembm.inc"
C
      REAL ALQ(2),BEQ(2),E,CS2THW,TNTHW,CTTHW,AE,BE,AM1,AM2,
     $EQ,ALR,Z,PHIZ,PROPZ,SIG,PCM,AMASS,ALL(2),BEL(2),
     $G,MSNE,TM2,TM3,TM4,TM5,TM6,AZJ,AZI,MEL,MER,
     $AEZS,BEZS,SR2,GP,AN,BN,AEZJS,BEZJS,SSXLAM,
     $TGG,TNN,TGN,TZN,AMWI,XS,YS,XC,YC,SINGL,SINGR,
     $COSGL,COSGR,XM,YM,THX,THY,XI,DEL,AMWISS(2),KK,
     $AMZIZ1,AMZIZ2,SIGLL,SIGRR,SIGLZ,SIGRZ,SSGT,SSGST,
     $FAC1,EZ0,BETA,EEL,EER,
     $FLEP,FLEM,FREP,FREM,SIGLR,SIGRL,PHIZLR,PHIZRL,
     $TM1LR,TM1RL,TZZRL,TZZLR,TGZLR,TGZRL,SIGZZL,SIGZZR,
     $FACLR,FACRL,RSH,JAC,ESTRUC,SH,SSFEL
      COMPLEX AEZ(4),BEZ(4),ZI,ZONE,WIJ
      INTEGER IS2UD(25),IUD(13),JS2JT(25),IQ1,IQ2,IFL1,IFL2,
     $IFLQ,IFM,I,IDQSS(25),MATCHL(18),IL2JS(18),IS2LN(18),
     $I1,I2,IL1,IL2,IDL1,IDL2,IZ,IZ1,IP,ITHZ(4),IDLSS(18),
     $IW2JS(4),IW1,JW1,JTW1,JTW2,IZ2JS(4),
     $IZ2,JTYPZ1,JTYPZ2
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
      DATA IDQSS/0,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2/
      DATA IDLSS/ISNEL,MSNEL,ISEL,MSEL,ISNML,MSNML,ISMUL,MSMUL,
     $ISNTL,MSNTL,ISTAU1,MSTAU1,ISER,MSER,ISMUR,MSMUR,
     $ISTAU2,MSTAU2/
      DATA IS2UD/0,1,1,2,2,2,2,1,1,2,2,1,1,1,1,2,2,2,2,1,1,2,2,1,1/
      DATA IUD/0,1,-1,2,-2,2,-2,1,-1,2,-2,1,-1/
      DATA JS2JT/1,
     $2,3,4,5,6,7,8,9,10,11,12,13,2,3,4,5,6,7,8,9,10,11,12,13/
      DATA MATCHL/2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/
      DATA IL2JS/34,35,36,37,38,39,40,41,42,43,44,45,46,47,
     $48,49,50,51/
      DATA IS2LN/1,1,2,2,1,1,2,2,1,1,2,2,2,2,2,2,2,2/
      DATA IW2JS/26,27,28,29/
      DATA IZ2JS/30,31,32,33/
      DATA ZONE,ZI/(1.,0.),(0.,1.)/
C
C          FUNCTIONS
      IF (IBREM) THEN
        SH=SHAT
        JAC=2*(1.-SHAT/SCM)*2*SQRT(SHAT)*(RSHMAX-RSHMIN)/SCM/(X1+X2)
      ELSE
        SH=SCM
      END IF
      PROPZ=(SH-AMZ**2)**2+AMZ**2*GAMZ**2
C
C          CONSTANTS
      RSH=SQRT(SH)
      EB=RSH/2.
      QSQBM=QSQ
      E=SQRT(4*PI*ALFAEM)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      BETA=ATAN(1./RV2V1)
      SR2=SQRT(2.)
      CS2THW=1.-SN2THW
      TNTHW=SQRT(SN2THW/CS2THW)
      CTTHW=1./TNTHW
      ALQ(1)=CTTHW/4.-5*TNTHW/12.
      BEQ(1)=-(CTTHW+TNTHW)/4.
      ALQ(2)=TNTHW/12.-CTTHW/4.
      BEQ(2)=-BEQ(1)
      ALL(1)=(CTTHW+TNTHW)/4.
      BEL(1)=-(CTTHW+TNTHW)/4.
      ALL(2)=(3*TNTHW-CTTHW)/4.
      BEL(2)=-BEL(1)
      AE=ALL(2)
      BE=BEL(2)
      AN=ALL(1)
      BN=BEL(1)
      FLEP=(1.+PLEP)/2.
      FLEM=(1.+PLEM)/2.
      FREP=(1.-PLEP)/2.
      FREM=(1.-PLEM)/2.
      MEL=AMASS(ISEL)
      MER=AMASS(ISER)
      MSNE=AMASS(ISNEL)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      AMWISS(1)=ABS(AMW1SS)
      AMWISS(2)=ABS(AMW2SS)
      DO 5 IZ=1,4
        ITHZ(IZ)=0
        IF (AMZISS(IZ).LT.0.) ITHZ(IZ)=1
        AEZ(IZ)=-1*ZI**(ITHZ(IZ)-1)*(-1)**(ITHZ(IZ)+1)*
     $      (G*ZMIXSS(3,IZ)+GP*ZMIXSS(4,IZ))/SR2
        BEZ(IZ)=-1*ZI**(ITHZ(IZ)-1)*SR2*GP*ZMIXSS(4,IZ)
5     CONTINUE
C
C          ENTRY
      SIG=0.
      SIGMA=0.
      NSIGS=0
      DO 10 I=1,MXSIGS
        SIGS(I)=0.
10    CONTINUE
C
C          First do squark pairs: IQ1 labels JETTYPE1.
C
      DO 100 IQ1=2,25
        IQ2=MATCH(IQ1,4)
        IF(.NOT.(GOQ(IQ1,1).AND.GOQ(IQ2,2))) GO TO 100
        IFL1=IDQSS(IQ1)
        IFL2=IDQSS(IQ2)
        AM1=AMASS(IFL1)
        AM2=AMASS(IFL2)
        IF((AM1+AM2).GE.RSH) GO TO 100
        IFLQ=IS2UD(IQ1)
        IF (IFLQ.EQ.1) THEN
          EQ=2./3.
        ELSE
          EQ=-1./3.
        END IF
C          Left squarks
        IF(IQ1.LE.9) THEN
          ALR=2*(ALQ(IFLQ)-BEQ(IFLQ))
C          Right squarks
        ELSEIF(IQ1.GE.14.AND.IQ1.LE.21) THEN
          ALR=2*(ALQ(IFLQ)+BEQ(IFLQ))
C          Mixed stops and sbottoms
        ELSEIF(IQ1.EQ.10.OR.IQ1.EQ.11) THEN
          ALR=2*(ALQ(IFLQ)-BEQ(IFLQ)*COS(2*THETAB))
        ELSEIF(IQ1.EQ.12.OR.IQ1.EQ.13) THEN
          ALR=2*(ALQ(IFLQ)-BEQ(IFLQ)*COS(2*THETAT))
        ELSEIF(IQ1.EQ.22.OR.IQ1.EQ.23) THEN
          ALR=2*(ALQ(IFLQ)+BEQ(IFLQ)*COS(2*THETAB))
        ELSEIF(IQ1.EQ.24.OR.IQ1.EQ.25) THEN
          ALR=2*(ALQ(IFLQ)+BEQ(IFLQ)*COS(2*THETAT))
        END IF
        PCM=.5*SQRT(SH-4.*AM1**2)
        IFM=ISIGN(1,IUD(JS2JT(IQ1)))
        IF (IFM.GT.0) THEN
          Z=CTH(1)
        ELSE
          Z=-CTH(1)
        END IF
C          Calculate d(sigma)/d(cos theta) in mb
        PHIZLR=2*E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE-BE)**2*
     $   SH-8*(AE-BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        PHIZRL=2*E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE+BE)**2*
     $   SH-8*(AE+BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        SIGLR=3*PCM**3/512./PI/EB**3*PHIZLR
        SIGRL=3*PCM**3/512./PI/EB**3*PHIZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,IQ1,IQ2)
100   CONTINUE
C        Mixed sbottom_1 and sbottom_2 production
      IF ((AMB1SS+AMB2SS).LT.RSH) THEN
        Z=CTH(1)
        PCM=SQRT(SSXLAM(SH,AMB1SS**2,AMB2SS**2))/2./RSH
        SIGLR=2*3*8*PI*ALFAEM**2*BEQ(2)**2*COS(THETAB)**2*
     $   SIN(THETAB)**2*(AE-BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIGRL=2*3*8*PI*ALFAEM**2*BEQ(2)**2*COS(THETAB)**2*
     $   SIN(THETAB)**2*(AE+BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(10,1).AND.GOQ(23,2)) THEN
          CALL SIGFIL(SIG,0,0,10,23)
        END IF
        IF(GOQ(23,1).AND.GOQ(10,2)) THEN
          CALL SIGFIL(SIG,0,0,23,10)
        END IF
        IF(GOQ(11,1).AND.GOQ(22,2)) THEN
          CALL SIGFIL(SIG,0,0,11,22)
        END IF
        IF(GOQ(22,1).AND.GOQ(11,2)) THEN
          CALL SIGFIL(SIG,0,0,22,11)
        END IF
      ENDIF
C        Mixed stop_1 and stop_2 production
      IF ((AMT1SS+AMT2SS).LT.RSH) THEN
        Z=CTH(1)
        PCM=SQRT(SSXLAM(SH,AMT1SS**2,AMT2SS**2))/2./RSH
        SIGLR=2*3*8*PI*ALFAEM**2*BEQ(1)**2*COS(THETAT)**2*
     $   SIN(THETAT)**2*(AE-BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIGRL=2*3*8*PI*ALFAEM**2*BEQ(1)**2*COS(THETAT)**2*
     $   SIN(THETAT)**2*(AE+BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(12,1).AND.GOQ(25,2)) THEN
          CALL SIGFIL(SIG,0,0,12,25)
        END IF
        IF(GOQ(25,1).AND.GOQ(12,2)) THEN
          CALL SIGFIL(SIG,0,0,25,12)
        END IF
        IF(GOQ(13,1).AND.GOQ(24,2)) THEN
          CALL SIGFIL(SIG,0,0,13,24)
        END IF
        IF(GOQ(24,1).AND.GOQ(13,2)) THEN
          CALL SIGFIL(SIG,0,0,24,13)
        END IF
      ENDIF
C
C          2nd and 3rd generation sleptons: IL1 labels JETTYPE1.
C
      DO 200 I=5,16
        I1=I
        IF (I1.GE.13) I1=I1+2
        I2=MATCHL(I1)
        IL1=IL2JS(I1)
        IL2=IL2JS(I2)
        IF(.NOT.(GOQ(IL1,1).AND.GOQ(IL2,2))) GO TO 200
        IDL1=IDLSS(I1)
        IDL2=IDLSS(I2)
        AM1=AMASS(IDL1)
        AM2=AMASS(IDL2)
        IF((AM1+AM2).GE.RSH) GO TO 200
        IFL1=IS2LN(I1)
        IFL2=IS2LN(I2)
        IF (IFL1.EQ.1) THEN
          EQ=0.
        ELSE
          EQ=-1.
        END IF
        IF (I1.EQ.15.OR.I1.EQ.16)  THEN
          ALR=2*(ALL(IFL1)+BEL(IFL1))
        ELSE IF (I1.GE.5.AND.I1.LE.10) THEN
          ALR=2*(ALL(IFL1)-BEL(IFL1))
        ELSE IF (I1.EQ.11.OR.I1.EQ.12) THEN
          ALR=2*(ALL(IFL1)-BEL(IFL1)*COS(2*THETAL))
        ELSE IF (I1.EQ.17.OR.I1.EQ.18) THEN
          ALR=2*(ALL(IFL1)+BEL(IFL1)*COS(2*THETAL))
        END IF
        PCM=.5*SQRT(SH-4.*AM1**2)
        IFM=ISIGN(1,IDL1)
        IF (IFM.GT.0) THEN
          Z=CTH(1)
         ELSE
          Z=-CTH(1)
        END IF
C          Calculate d(sigma)/d(cos theta) in mb
        PHIZLR=2*E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE-BE)**2*
     $   SH-8*(AE-BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        PHIZRL=2*E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE+BE)**2*
     $   SH-8*(AE+BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        SIGLR=PCM**3/512./PI/EB**3*PHIZLR
        SIGRL=PCM**3/512./PI/EB**3*PHIZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,IL1,IL2)
200   CONTINUE
C        Mixed stau_1 and stau_2 production
      IF ((AML1SS+AML2SS).LT.RSH) THEN
        Z=CTH(1)
        PCM=SQRT(SSXLAM(SH,AML1SS**2,AML2SS**2))/2./RSH
        SIGLR=2*8*PI*ALFAEM**2*BEL(2)**2*COS(THETAL)**2*
     $   SIN(THETAL)**2*(AE-BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIGRL=2*8*PI*ALFAEM**2*BEL(2)**2*COS(THETAL)**2*
     $   SIN(THETAL)**2*(AE+BE)**2*PCM**3*(1.-Z**2)/RSH/PROPZ
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(44,1).AND.GOQ(51,2)) THEN
          CALL SIGFIL(SIG,0,0,44,51)
        END IF
        IF(GOQ(51,1).AND.GOQ(44,2)) THEN
          CALL SIGFIL(SIG,0,0,51,44)
        END IF
        IF(GOQ(45,1).AND.GOQ(50,2)) THEN
          CALL SIGFIL(SIG,0,0,45,50)
        END IF
        IF(GOQ(50,1).AND.GOQ(45,2)) THEN
          CALL SIGFIL(SIG,0,0,50,45)
        END IF
      ENDIF
C
C         Next do 1st generation sleptons
C
C         Sneutrino_e pairs
      DO 210 I1=1,2
        I2=MATCHL(I1)
        IL1=IL2JS(I1)
        IL2=IL2JS(I2)
        IF(.NOT.(GOQ(IL1,1).AND.GOQ(IL2,2))) GO TO 210
        MSNE=AMASS(ISNEL)
        IF((2*MSNE).GE.RSH) GO TO 210
        IF (I1.EQ.1) THEN
          Z=CTH(1)
        ELSE
          Z=-CTH(1)
        END IF
        PCM=.5*SQRT(SH-4*MSNE**2)
        TM1LR=32*E**4*(AN-BN)**2*(AE-BE)**2/PROPZ
        TM1RL=32*E**4*(AN-BN)**2*(AE+BE)**2/PROPZ
        TM2=8*G**4*SIN(GAMMAR)**4/(2*EB*(EB-PCM*Z)+AMW1SS**2-MSNE**2)**2
        TM3=8*G**4*COS(GAMMAR)**4/(2*EB*(EB-PCM*Z)+AMW2SS**2-MSNE**2)**2
        TM4=-32*E**2*(AN-BN)*G**2*SIN(GAMMAR)**2*(SH-AMZ**2)*(AE-BE)/
     $  PROPZ/(2*EB*(EB-PCM*Z)+AMW1SS**2-MSNE**2)
        TM5=-32*E**2*(AN-BN)*G**2*COS(GAMMAR)**2*(SH-AMZ**2)*(AE-BE)/
     $  PROPZ/(2*EB*(EB-PCM*Z)+AMW2SS**2-MSNE**2)
        TM6=16*G**4*SIN(GAMMAR)**2*COS(GAMMAR)**2/
     $  (2*EB*(EB-PCM*Z)+AMW1SS**2-MSNE**2)/
     $  (2*EB*(EB-PCM*Z)+AMW2SS**2-MSNE**2)
        SIGLR=2*PCM**3*EB*(1.-Z**2)/128./PI/SH*
     $   (TM1LR+TM2+TM3+TM4+TM5+TM6)
        SIGRL=2*PCM**3*EB*(1.-Z**2)/128./PI/SH*TM1RL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,IL1,IL2)
210   CONTINUE
C         E_L~ pairs
      DO 220 I1=3,4
        I2=MATCHL(I1)
        IL1=IL2JS(I1)
        IL2=IL2JS(I2)
        IF(.NOT.(GOQ(IL1,1).AND.GOQ(IL2,2))) GO TO 220
        IF(2*MEL.GE.RSH) GO TO 220
        PCM=.5*SQRT(SH-4.*MEL**2)
        EQ=-1.
        ALR=2*(AE-BE)
        IF (I1.EQ.3) THEN
          Z=CTH(1)
        ELSE
          Z=-CTH(1)
        END IF
        PHIZLR=E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE-BE)**2*
     $   SH-8*(AE-BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        PHIZRL=E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE+BE)**2*
     $   SH-8*(AE+BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        DO 221 IZ1=1,4
          AEZS=AEZ(IZ1)*CONJG(AEZ(IZ1))
          PHIZLR=PHIZLR+2*AEZS**2*SH*(1.-Z**2)/(2*EB*(EB-PCM*Z)-
     $    MEL**2+AMZISS(IZ1)**2)**2-4*E**2*(1.-Z**2)*AEZS/
     $    (2*EB*(EB-PCM*Z)-MEL**2+AMZISS(IZ1)**2)*(2.+(AE-BE)*ALR*
     $    SH*(SH-AMZ**2)/PROPZ)
          IF (IZ1.LE.3) THEN
            DO 222 IP=IZ1+1,4
              AEZJS=AEZ(IP)*CONJG(AEZ(IP))
              PHIZLR=PHIZLR+4*AEZS*AEZJS*SH*(1.-Z**2)/
     $        (2*EB*(EB-PCM*Z)-MEL**2+AMZISS(IZ1)**2)/
     $        (2*EB*(EB-PCM*Z)-MEL**2+AMZISS(IP)**2)
222         CONTINUE
          END IF
221     CONTINUE
        SIGLR=2*PCM**3/512./PI/EB**3*PHIZLR
        SIGRL=2*PCM**3/512./PI/EB**3*PHIZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,IL1,IL2)
220   CONTINUE
C         E_R~ pairs
      DO 230 I1=13,14
        I2=MATCHL(I1)
        IL1=IL2JS(I1)
        IL2=IL2JS(I2)
        IF(.NOT.(GOQ(IL1,1).AND.GOQ(IL2,2))) GO TO 230
        IF(2*MER.GE.RSH) GO TO 230
        PCM=.5*SQRT(SH-4.*MER**2)
        EQ=-1.
        ALR=2*(AE+BE)
        IF (I1.EQ.13) THEN
          Z=CTH(1)
        ELSE
          Z=-CTH(1)
        END IF
        PHIZLR=E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE-BE)**2*
     $   SH-8*(AE-BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        PHIZRL=E**4*(1.-Z**2)*(8*EQ**2/SH+(2*ALR**2*(AE+BE)**2*
     $   SH-8*(AE+BE)*EQ*ALR*(SH-AMZ**2))/PROPZ)
        DO 231 IZ1=1,4
          BEZS=BEZ(IZ1)*CONJG(BEZ(IZ1))
          PHIZRL=PHIZRL+2*BEZS**2*SH*(1.-Z**2)/(2*EB*(EB-PCM*Z)-
     $    MER**2+AMZISS(IZ1)**2)**2-4*E**2*(1.-Z**2)*BEZS/
     $    (2*EB*(EB-PCM*Z)-MER**2+AMZISS(IZ1)**2)*(2.+(AE+BE)*ALR*
     $    SH*(SH-AMZ**2)/PROPZ)
          IF (IZ1.LE.3) THEN
            DO 232 IP=IZ1+1,4
              BEZJS=BEZ(IP)*CONJG(BEZ(IP))
              PHIZRL=PHIZRL+4*BEZS*BEZJS*SH*(1.-Z**2)/
     $        (2*EB*(EB-PCM*Z)-MER**2+AMZISS(IZ1)**2)/
     $        (2*EB*(EB-PCM*Z)-MER**2+AMZISS(IP)**2)
232         CONTINUE
          END IF
231     CONTINUE
        SIGLR=2*PCM**3/512./PI/EB**3*PHIZLR
        SIGRL=2*PCM**3/512./PI/EB**3*PHIZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,IL1,IL2)
230   CONTINUE
C         E_L~+E_R~bar and E_R~+E_L~bar pairs; now has MEL =/ MER !
      IF((MEL+MER).GE.RSH) GO TO 270
      IF(GOQ(36,1).AND.GOQ(47,2)) THEN
        PCM=SQRT(SSXLAM(SH,MEL**2,MER**2))/4./EB
        EEL=SQRT(PCM**2+MEL**2)
        Z=CTH(1)
        PHIZ=0.
        DO 241 IZ1=1,4
          BEZS=BEZ(IZ1)*CONJG(BEZ(IZ1))
          AEZS=AEZ(IZ1)*CONJG(AEZ(IZ1))
          AZI=(AMZISS(IZ1)**2-MEL**2)/2./EB
          PHIZ=PHIZ+AEZS*BEZS*AMZISS(IZ1)**2/(EEL-PCM*Z+AZI)**2
          IF (IZ1.LE.3) THEN
            DO 242 IP=IZ1+1,4
              AZJ=(AMZISS(IP)**2-MEL**2)/2./EB
              PHIZ=PHIZ+2*ABS(AMZISS(IZ1)*AMZISS(IP))*
     $        REAL(AEZ(IZ1)*CONJG(AEZ(IP))*CONJG(BEZ(IZ1))*BEZ(IP))/
     $        (EEL-PCM*Z+AZI)/(EEL-PCM*Z+AZJ)
242         CONTINUE
          END IF
241     CONTINUE
        SIG=4*PCM/128./PI/SH/EB*PHIZ
        SIG=FLEM*FLEP*SIG*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,36,47)
      ENDIF
      IF(GOQ(46,1).AND.GOQ(37,2)) THEN
        PCM=SQRT(SSXLAM(SH,MEL**2,MER**2))/4./EB
        EER=SQRT(PCM**2+MER**2)
        Z=CTH(1)
        PHIZ=0.
        DO 243 IZ1=1,4
          BEZS=BEZ(IZ1)*CONJG(BEZ(IZ1))
          AEZS=AEZ(IZ1)*CONJG(AEZ(IZ1))
          AZI=(AMZISS(IZ1)**2-MER**2)/2./EB
          PHIZ=PHIZ+AEZS*BEZS*AMZISS(IZ1)**2/(EER-PCM*Z+AZI)**2
          IF (IZ1.LE.3) THEN
            DO 244 IP=IZ1+1,4
              AZJ=(AMZISS(IP)**2-MER**2)/2./EB
              PHIZ=PHIZ+2*ABS(AMZISS(IZ1)*AMZISS(IP))*
     $        REAL(AEZ(IZ1)*CONJG(AEZ(IP))*CONJG(BEZ(IZ1))*BEZ(IP))/
     $        (EER-PCM*Z+AZI)/(EER-PCM*Z+AZJ)
244         CONTINUE
          END IF
243     CONTINUE
        SIG=4*PCM/128./PI/SH/EB*PHIZ
        SIG=FREM*FREP*SIG*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,46,37)
      ENDIF
C         E_R~bar+E_L~ and E_L~bar+E_R~ pairs; now assumes MEL =/ MER !
      IF(GOQ(47,1).AND.GOQ(36,2)) THEN
        PCM=SQRT(SSXLAM(SH,MEL**2,MER**2))/4./EB
        EEL=SQRT(PCM**2+MEL**2)
        Z=-CTH(1)
        PHIZ=0.
        DO 251 IZ1=1,4
          BEZS=BEZ(IZ1)*CONJG(BEZ(IZ1))
          AEZS=AEZ(IZ1)*CONJG(AEZ(IZ1))
          AZI=(AMZISS(IZ1)**2-MEL**2)/2./EB
          PHIZ=PHIZ+AEZS*BEZS*AMZISS(IZ1)**2/(EEL-PCM*Z+AZI)**2
          IF (IZ1.LE.3) THEN
            DO 252 IP=IZ1+1,4
              AZJ=(AMZISS(IP)**2-MEL**2)/2./EB
              PHIZ=PHIZ+2*ABS(AMZISS(IZ1)*AMZISS(IP))*
     $        REAL(AEZ(IZ1)*CONJG(AEZ(IP))*CONJG(BEZ(IZ1))*BEZ(IP))/
     $        (EEL-PCM*Z+AZI)/(EEL-PCM*Z+AZJ)
252         CONTINUE
          END IF
251     CONTINUE
        SIG=4*PCM/128./PI/SH/EB*PHIZ
        SIG=FLEM*FLEP*SIG*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,47,36)
      ENDIF
      IF(GOQ(37,1).AND.GOQ(46,2)) THEN
        PCM=SQRT(SSXLAM(SH,MEL**2,MER**2))/4./EB
        EER=SQRT(PCM**2+MER**2)
        Z=-CTH(1)
        PHIZ=0.
        DO 253 IZ1=1,4
          BEZS=BEZ(IZ1)*CONJG(BEZ(IZ1))
          AEZS=AEZ(IZ1)*CONJG(AEZ(IZ1))
          AZI=(AMZISS(IZ1)**2-MER**2)/2./EB
          PHIZ=PHIZ+AEZS*BEZS*AMZISS(IZ1)**2/(EER-PCM*Z+AZI)**2
          IF (IZ1.LE.3) THEN
            DO 254 IP=IZ1+1,4
              AZJ=(AMZISS(IP)**2-MER**2)/2./EB
              PHIZ=PHIZ+2*ABS(AMZISS(IZ1)*AMZISS(IP))*
     $        REAL(AEZ(IZ1)*CONJG(AEZ(IP))*CONJG(BEZ(IZ1))*BEZ(IP))/
     $        (EER-PCM*Z+AZI)/(EER-PCM*Z+AZJ)
254         CONTINUE
          END IF
253     CONTINUE
        SIG=4*PCM/128./PI/SH/EB*PHIZ
        SIG=FREM*FREP*SIG*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,37,46)
      ENDIF
270   CONTINUE
C
C          Chargino pair production
C
      DO 300 IW1=1,4
        JW1=(IW1+1)/2
        AMWI=ABS(AMWISS(JW1))
        JTW1=IW2JS(IW1)
        JTW2=IW2JS(MATCHL(IW1))
        IF (.NOT.(GOQ(JTW1,1).AND.GOQ(JTW2,2))) GO TO 300
        IF((2*AMWI).GE.RSH) GO TO 300
        PCM=SQRT(SSXLAM(SH,AMWI**2,AMWI**2))/4./EB
        Z=CTH(1)
        IF (IW1.EQ.1.OR.IW1.EQ.3) Z=-CTH(1)
        SINGR=SIN(GAMMAR)
        COSGR=COS(GAMMAR)
        SINGL=SIN(GAMMAL)
        COSGL=COS(GAMMAL)
        XC=1.-(COSGL**2+COSGR**2)/4./CS2THW
        YC=(COSGR**2-COSGL**2)/4./CS2THW
        XS=1.-(SINGL**2+SINGR**2)/4./CS2THW
        YS=(SINGR**2-SINGL**2)/4./CS2THW
        IF (IW1.GE.3) THEN
          XC=XS
          YC=YS
          SINGR=COSGR
        END IF
        TGG=16*E**4/SH*(EB**2*(1.+Z**2)+AMWI**2*(1.-Z**2))
        TZZLR=16*E**4*CTTHW**2*SH/PROPZ*((XC**2+YC**2)*(AE-BE)**2*
     $  (EB**2*(1.+Z**2)+AMWI**2*(1.-Z**2))-
     $  2*YC**2*(AE-BE)**2*AMWI**2+4*XC*YC*(AE-BE)**2*EB*PCM*Z)
        TZZRL=16*E**4*CTTHW**2*SH/PROPZ*((XC**2+YC**2)*(AE+BE)**2*
     $  (EB**2*(1.+Z**2)+AMWI**2*(1.-Z**2))-
     $  2*YC**2*(AE+BE)**2*AMWI**2-4*XC*YC*(AE+BE)**2*EB*PCM*Z)
        TGZLR=-32*E**4*CTTHW*(SH-AMZ**2)/PROPZ*((AE-BE)*XC*
     $  (EB**2*(1.+Z**2)+AMWI**2*(1.-Z**2))-2*(BE-AE)*YC*EB*PCM*Z)
        TGZRL=-32*E**4*CTTHW*(SH-AMZ**2)/PROPZ*((AE+BE)*XC*
     $  (EB**2*(1.+Z**2)+AMWI**2*(1.-Z**2))-2*(BE+AE)*YC*EB*PCM*Z)
        TNN=2*E**4*SINGR**4*SH*(EB-PCM*Z)**2/SN2THW**2/
     $  (EB**2+PCM**2-2*EB*PCM*Z+MSNE**2)**2
        TGN=-8*E**4*SINGR**2*((EB-PCM*Z)**2+AMWI**2)/SN2THW/
     $  (EB**2+PCM**2-2*EB*PCM*Z+MSNE**2)
        TZN=8*E**4*CTTHW*SINGR**2*(SH-AMZ**2)*(AE-BE)*SH/
     $  SN2THW/PROPZ*((XC-YC)*((EB-PCM*Z)**2+AMWI**2)+2*YC*AMWI**2)/
     $  (EB**2+PCM**2-2*EB*PCM*Z+MSNE**2)
        SIGLR=2*PCM/128./PI/SH/EB*(TGG+TZZLR+TGZLR+TNN+TGN+TZN)
        SIGRL=2*PCM/128./PI/SH/EB*(TGG+TZZRL+TGZRL)
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,JTW1,JTW2)
300   CONTINUE
C
C     Chargino_1 + chargino_2 pair production
      IF((ABS(AMW1SS)+ABS(AMW2SS)).GE.RSH) GO TO 340
      PCM=SQRT(SSXLAM(SH,AMW1SS**2,AMW2SS**2))/4./EB
      XC=(THX*SIN(GAMMAL)*COS(GAMMAL)-THY*SIN(GAMMAR)*COS(GAMMAR))/2.
      YC=(THX*SIN(GAMMAL)*COS(GAMMAL)+THY*SIN(GAMMAR)*COS(GAMMAR))/2.
      DEL=(AMW2SS**2-AMW1SS**2)/4./EB
      XI=-1.*SIGN(1.,AMWISS(1))*SIGN(1.,AMWISS(2))
      IF (.NOT.(GOQ(27,1).AND.GOQ(28,2))) GO TO 310
        Z=CTH(1)
        TZZLR=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE-BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE-BE)**2*ABS(AMW1SS*AMW2SS)+
     $  4*XC*YC*(AE-BE)**2*EB*PCM*Z)
        TZZRL=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE+BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE+BE)**2*ABS(AMW1SS*AMW2SS)-
     $  4*XC*YC*(AE+BE)**2*EB*PCM*Z)
        TNN=2*SIN(GAMMAR)**2*COS(GAMMAR)**2*((EB-PCM*Z)**2-DEL**2)/
     $  SN2THW**2/(2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)**2
        TZN=-4*THY*(CTTHW+TNTHW)*SIN(GAMMAR)*COS(GAMMAR)*(SH-AMZ**2)
     $  *(AE-BE)/SN2THW/PROPZ*((XC-YC)*((EB-PCM*Z)**2-DEL**2-
     $  XI*ABS(AMW1SS*AMW2SS))+2*XC*XI*ABS(AMW1SS*AMW2SS))/
     $  (2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)
        SIGLR=2*E**4*PCM/128./PI/EB*(TZZLR+TNN+TZN)
        SIGRL=2*E**4*PCM/128./PI/EB*TZZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,27,28)
310   CONTINUE
      IF (.NOT.(GOQ(28,1).AND.GOQ(27,2))) GO TO 320
        Z=-CTH(1)
        TZZLR=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE-BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE-BE)**2*ABS(AMW1SS*AMW2SS)+
     $  4*XC*YC*(AE-BE)**2*EB*PCM*Z)
        TZZRL=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE+BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE+BE)**2*ABS(AMW1SS*AMW2SS)-
     $  4*XC*YC*(AE+BE)**2*EB*PCM*Z)
        TNN=2*SIN(GAMMAR)**2*COS(GAMMAR)**2*((EB-PCM*Z)**2-DEL**2)/
     $  SN2THW**2/(2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)**2
        TZN=-4*THY*(CTTHW+TNTHW)*SIN(GAMMAR)*COS(GAMMAR)*(SH-AMZ**2)
     $  *(AE-BE)/SN2THW/PROPZ*((XC-YC)*((EB-PCM*Z)**2-DEL**2-
     $  XI*ABS(AMW1SS*AMW2SS))+2*XC*XI*ABS(AMW1SS*AMW2SS))/
     $  (2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)
        SIGLR=2*E**4*PCM/128./PI/EB*(TZZLR+TNN+TZN)
        SIGRL=2*E**4*PCM/128./PI/EB*TZZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,28,27)
320   CONTINUE
      IF (.NOT.(GOQ(29,1).AND.GOQ(26,2))) GO TO 330
        Z=CTH(1)
        TZZLR=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE-BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE-BE)**2*ABS(AMW1SS*AMW2SS)+
     $  4*XC*YC*(AE-BE)**2*EB*PCM*Z)
        TZZRL=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE+BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE+BE)**2*ABS(AMW1SS*AMW2SS)-
     $  4*XC*YC*(AE+BE)**2*EB*PCM*Z)
        TNN=2*SIN(GAMMAR)**2*COS(GAMMAR)**2*((EB-PCM*Z)**2-DEL**2)/
     $  SN2THW**2/(2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)**2
        TZN=-4*THY*(CTTHW+TNTHW)*SIN(GAMMAR)*COS(GAMMAR)*(SH-AMZ**2)
     $  *(AE-BE)/SN2THW/PROPZ*((XC-YC)*((EB-PCM*Z)**2-DEL**2-
     $  XI*ABS(AMW1SS*AMW2SS))+2*XC*XI*ABS(AMW1SS*AMW2SS))/
     $  (2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)
        SIGLR=2*E**4*PCM/128./PI/EB*(TZZLR+TNN+TZN)
        SIGRL=2*E**4*PCM/128./PI/EB*TZZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,29,26)
330   CONTINUE
      IF (.NOT.(GOQ(26,1).AND.GOQ(29,2))) GO TO 340
        Z=-CTH(1)
        TZZLR=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE-BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE-BE)**2*ABS(AMW1SS*AMW2SS)+
     $  4*XC*YC*(AE-BE)**2*EB*PCM*Z)
        TZZRL=4*(CTTHW+TNTHW)**2/PROPZ*((XC**2+YC**2)*(AE+BE)**2*
     $  (EB**2+PCM**2*Z**2-DEL**2-XI*ABS(AMW1SS*AMW2SS))+
     $  2*XC**2*XI*(AE+BE)**2*ABS(AMW1SS*AMW2SS)-
     $  4*XC*YC*(AE+BE)**2*EB*PCM*Z)
        TNN=2*SIN(GAMMAR)**2*COS(GAMMAR)**2*((EB-PCM*Z)**2-DEL**2)/
     $  SN2THW**2/(2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)**2
        TZN=-4*THY*(CTTHW+TNTHW)*SIN(GAMMAR)*COS(GAMMAR)*(SH-AMZ**2)
     $  *(AE-BE)/SN2THW/PROPZ*((XC-YC)*((EB-PCM*Z)**2-DEL**2-
     $  XI*ABS(AMW1SS*AMW2SS))+2*XC*XI*ABS(AMW1SS*AMW2SS))/
     $  (2*EB*(EB-DEL)-2*EB*PCM*Z+MSNE**2-AMW1SS**2)
        SIGLR=2*E**4*PCM/128./PI/EB*(TZZLR+TNN+TZN)
        SIGRL=2*E**4*PCM/128./PI/EB*TZZRL
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        CALL SIGFIL(SIG,0,0,26,29)
340   CONTINUE
C
C         Neutralino pair production
C
      DO 400 IZ1=1,4
        AMZIZ1=ABS(AMZISS(IZ1))
        JTYPZ1=IZ2JS(IZ1)
        DO 410 IZ2=1,4
          AMZIZ2=ABS(AMZISS(IZ2))
          JTYPZ2=IZ2JS(IZ2)
          IF(.NOT.(GOQ(JTYPZ1,1).AND.GOQ(JTYPZ2,2))) GO TO 410
          IF((AMZIZ1+AMZIZ2).GE.RSH) GO TO 410
          WIJ=SQRT(G**2+GP**2)*ZI**(ITHZ(IZ2))*(-ZI)**(ITHZ(IZ1))*
     $    (ZMIXSS(1,IZ1)*ZMIXSS(1,IZ2)-ZMIXSS(2,IZ1)*
     $    ZMIXSS(2,IZ2))/4.
          KK=SQRT(SH*SH+(AMZIZ1**2-AMZIZ2**2)**2-2*SH*
     $    (AMZIZ1**2+AMZIZ2**2))/4./EB
          Z=CTH(1)
          SIGLL=2*AEZ(IZ1)*CONJG(AEZ(IZ1))*AEZ(IZ2)*CONJG(AEZ(IZ2))*
     $    SSGT(SH,MEL,Z,IZ1,IZ2)
          SIGRR=2*BEZ(IZ1)*CONJG(BEZ(IZ1))*BEZ(IZ2)*CONJG(BEZ(IZ2))*
     $    SSGT(SH,MER,Z,IZ1,IZ2)
          SIGZZL=4*E**2*WIJ*CONJG(WIJ)*(AE-BE)**2*
     $    (SH*SH-(AMZIZ1**2-AMZIZ2**2)**2+4*(-1.)**(ITHZ(IZ1)+
     $    ITHZ(IZ2)+1)*SH*AMZIZ1*AMZIZ2+4*SH*KK*KK*Z*Z)/PROPZ
          SIGZZR=4*E**2*WIJ*CONJG(WIJ)*(AE+BE)**2*
     $    (SH*SH-(AMZIZ1**2-AMZIZ2**2)**2+4*(-1.)**(ITHZ(IZ1)+
     $    ITHZ(IZ2)+1)*SH*AMZIZ1*AMZIZ2+4*SH*KK*KK*Z*Z)/PROPZ
          SIGLZ=-E*(AE-BE)*(SH-AMZ**2)/2./PROPZ*
     $    (REAL(WIJ*CONJG(AEZ(IZ1))*AEZ(IZ2))*
     $    SSGST(SH,MEL,Z,IZ1,IZ2)+(-1.)**(ITHZ(IZ1)+ITHZ(IZ2))*
     $    REAL(WIJ*AEZ(IZ1)*CONJG(AEZ(IZ2)))*
     $    SSGST(SH,MEL,-Z,IZ1,IZ2))
          SIGRZ=-E*(-1.)**(ITHZ(IZ1)+ITHZ(IZ2)+1)*
     $    (AE+BE)*(SH-AMZ**2)/2./PROPZ*
     $    (REAL(WIJ*CONJG(BEZ(IZ1))*BEZ(IZ2))*
     $    SSGST(SH,MER,Z,IZ1,IZ2)+(-1.)**(ITHZ(IZ1)+ITHZ(IZ2))*
     $    REAL(WIJ*BEZ(IZ1)*CONJG(BEZ(IZ2)))*
     $    SSGST(SH,MER,-Z,IZ1,IZ2))
          SIGLR=2*KK/16./PI/SH/SQRT(SH)*(SIGLL+SIGZZL+SIGLZ)
          SIGRL=2*KK/16./PI/SH/SQRT(SH)*(SIGRR+SIGZZR+SIGRZ)
C         BELOW FACTOR OF 2 FOR ID PARTICLES AND JETTYP SWITCH
          SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
          IF (IBREM.AND..NOT.IBEAM) THEN
            SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
          ELSE IF (IBEAM) THEN
            SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
          END IF
          CALL SIGFIL(SIG,0,0,JTYPZ1,JTYPZ2)
410     CONTINUE
400   CONTINUE
C
C     Higgs boson mechanisms
C
C          E+ E- --> Z H_L; symmetric in cos(theta)
      IF((AMZ+AMHL).LT.RSH) THEN
        FACLR=E**2*G**2*(SIN(ALFAH+BETA))**2*(AE-BE)**2/CS2THW
        FACRL=E**2*G**2*(SIN(ALFAH+BETA))**2*(AE+BE)**2/CS2THW
        Z=CTH(1)
        PCM=SQRT(SSXLAM(SH,AMZ**2,AMHL**2))/4./EB
        EZ0=SQRT(PCM**2+AMZ**2)
        FAC1=AMZ**2+EZ0**2-PCM**2*Z**2
        SIGLR=2*FACLR/32./PI/PROPZ/SQRT(SH)*PCM*FAC1
        SIGRL=2*FACRL/32./PI/PROPZ/SQRT(SH)*PCM*FAC1
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(80,1).AND.GOQ(81,2)) CALL SIGFIL(SIG,0,0,80,81)
        IF(GOQ(81,1).AND.GOQ(80,2)) CALL SIGFIL(SIG,0,0,81,80)
      ENDIF
C          E+ E- --> Z H_H; symmetric in cos(theta)
      IF((AMZ+AMHH).LT.RSH) THEN
        FACLR=E**2*G**2*(COS(ALFAH+BETA))**2*(AE-BE)**2/CS2THW
        FACRL=E**2*G**2*(COS(ALFAH+BETA))**2*(AE+BE)**2/CS2THW
        Z=CTH(1)
        PCM=SQRT(SSXLAM(SH,AMZ**2,AMHH**2))/4./EB
        EZ0=SQRT(PCM**2+AMZ**2)
        FAC1=AMZ**2+EZ0**2-PCM**2*Z**2
        SIGLR=2*FACLR/32./PI/PROPZ/SQRT(SH)*PCM*FAC1
        SIGRL=2*FACRL/32./PI/PROPZ/SQRT(SH)*PCM*FAC1
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(80,1).AND.GOQ(82,2)) CALL SIGFIL(SIG,0,0,80,82)
        IF(GOQ(82,1).AND.GOQ(80,2)) CALL SIGFIL(SIG,0,0,82,80)
      ENDIF
C          E+ E- --> H_P H_L; symmetric in cos(theta)
      IF((AMHA+AMHL).LT.RSH) THEN
        PCM=SQRT(SSXLAM(SH,AMHA**2,AMHL**2))/4./EB
        Z=CTH(1)
        FAC1=PCM**3*(1.-Z**2)
        FACLR=E**4*(COS(ALFAH+BETA))**2*(AE-BE)**2*FAC1
        FACRL=E**4*(COS(ALFAH+BETA))**2*(AE+BE)**2*FAC1
        SIGLR=2*FACLR/32./PI/SQRT(SH)/SN2THW/CS2THW/PROPZ
        SIGRL=2*FACRL/32./PI/SQRT(SH)/SN2THW/CS2THW/PROPZ
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(81,1).AND.GOQ(83,2)) CALL SIGFIL(SIG,0,0,81,83)
        IF(GOQ(83,1).AND.GOQ(81,2)) CALL SIGFIL(SIG,0,0,83,81)
      ENDIF
C          E+ E- --> H_P H_H; SYMMETRIC IN COS(THETA)
      IF((AMHA+AMHH).LT.RSH) THEN
        PCM=SQRT(SSXLAM(SH,AMHA**2,AMHH**2))/4./EB
        Z=CTH(1)
        FAC1=PCM**3*(1.-Z**2)
        FACLR=E**4*(SIN(ALFAH+BETA))**2*(AE-BE)**2*FAC1
        FACRL=E**4*(SIN(ALFAH+BETA))**2*(AE+BE)**2*FAC1
        SIGLR=2*FACLR/32./PI/SQRT(SH)/SN2THW/CS2THW/PROPZ
        SIGRL=2*FACRL/32./PI/SQRT(SH)/SN2THW/CS2THW/PROPZ
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(82,1).AND.GOQ(83,2)) CALL SIGFIL(SIG,0,0,82,83)
        IF(GOQ(83,1).AND.GOQ(82,2)) CALL SIGFIL(SIG,0,0,83,82)
      ENDIF
C          E+ E- --> H^+ H^-; symmetric in cos(theta)
      IF((2*AMHC).LT.RSH) THEN
        PCM=SQRT(SSXLAM(SH,AMHC**2,AMHC**2))/4./EB
        Z=CTH(1)
        FAC1=PCM**3*(1.-Z**2)
        FACLR=FAC1*(1./SH**2+(2*SN2THW-1.)**2/SN2THW/CS2THW*
     $(AE-BE)**2/4./PROPZ+(2*SN2THW-1.)*(AE-BE)*(SH-AMZ**2)/SH/
     $SQRT(SN2THW*CS2THW)/PROPZ)
        FACRL=FAC1*(1./SH**2+(2*SN2THW-1.)**2/SN2THW/CS2THW*
     $(AE+BE)**2/4./PROPZ+(2*SN2THW-1.)*(AE+BE)*(SH-AMZ**2)/SH/
     $SQRT(SN2THW*CS2THW)/PROPZ)
        SIGLR=2*E**4*FACLR/8./PI/SQRT(SH)
        SIGRL=2*E**4*FACRL/8./PI/SQRT(SH)
        SIG=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS/2.
        IF (IBREM.AND..NOT.IBEAM) THEN
          SIG=SIG*ESTRUC(X1,QSQ)*ESTRUC(X2,QSQ)*JAC
        ELSE IF (IBEAM) THEN
          SIG=SIG*SSFEL(X1,0)*SSFEL(X2,0)*JAC
        END IF
        IF(GOQ(84,1).AND.GOQ(85,2)) CALL SIGFIL(SIG,0,0,84,85)
        IF(GOQ(85,1).AND.GOQ(84,2)) CALL SIGFIL(SIG,0,0,85,84)
      ENDIF
C
      RETURN
      END
