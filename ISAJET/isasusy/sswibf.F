#include "isajet/pilot.h"
        SUBROUTINE SSWIBF
C-----------------------------------------------------------------------
C       This subroutine calculates the chargino (wi) subset of
C       SSWZBF, which was too long.
C       Valid for all scalar masses (functions in double precision)
C       Includes Higgs sector radiative corrections (Aug. 31)
C
C       Auxiliary functions are called SSWxyi, SSZxyi, where normally
C       x indicates the SUSY particle, y the SM particle(s), and i is
C       a counter.
C
C       Part of Baer's GAUGBF
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
#include "isajet/sstmp.inc"
#include "isajet/sspols.inc"
C
      EXTERNAL SSZWF1,SSZZF1,SSZZF2,SSZZF3,SSWZF1,SSWZF2,SSWZF3
     $,SSWZF4,SSWZF5,SSWZF6,SSWZF7,SSWWF1,SSZZF4,SSZZF5,SSGX1
     $,SSGX2,SSGX8
C
      REAL AUI(4),BUI(4),ADI(4),BDI(4),ALI(4),BLI(4),ANI(4),BNI(4)
     $,WIJ(4,4),AUWI(4),ADWI(4),ANWI(4),ALWI(4),XIM(4),YIM(4)
     $,XIP(4),YIP(4),SNIJ(4,4),XLIJ(4,4),HIJ(4,4)
     $,V1I(4),V2I(4),V3I(4),V4I(4),XHIJ(4,4),XPIJ(4,4),AMWISS(2)
C
      INTEGER ISZ(4),THJZ
C
      REAL MWIW,SL,PP,SP,PL,PH,SH,COSB,SINB,COSA,SINA
     $,FACTOR,MZJZ,ULIM,XINTGL,COSBE,SINBE
     $,UPPER,CONST
      REAL W21ZU,W21ZN,W21ZL,W21ZD,W21U,W21D,W21S,W21C,W21N1
     $,W21N2,W21N3,W21E,W21M,W21L,STHW,CTHW
      REAL T3,XI2,CC,PSIINT,T2,T1,PHIINT,XI1,EF,A,Z,B
     $,TANB,FB,FT,SR2,G,PI,GP,FL
     $,MW2,SNW1,MW1,YM,BE,SNW2,XM,THX,THY
     $,BTN,APD,APL,APU,BTD,BTL,APN,BTU,Y,MZ1,FPI
      REAL TANW,COTW,XWINO,YWINO,SNIW,SNJZ
      REAL SSXINT,SSXLAM
      REAL WID,TERM1,TERM2,TERM3,TERM4,E,TERMH
      REAL FACT,ALJZ1,ALJZ2,BEJZ1,BEJZ2
      REAL ALIW1,ALIW2,AHCJZ,BHCJZ,TERMW,TERMN,TERM12,TERMN1,
     $TERMN2,TERMWN,TERMW1,TERMW2,TERMH1,TERMH2,TERMHN
      REAL XIPM,YIPM,COSL,SINL,BPWI(2),BPLWI(2)
      REAL BWI(2),AS,COST,SINT
      REAL POLNL,POLNR,POL1L,POL1R,POL2L,POL2R,POL12L,POL12R,
     $POLN1L,POLN1R,POLN2L,POLN2R
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      REAL SUALFE,MTAMTA,MTAMB,MTAMZ,AMPL,AMPI
      REAL FUDGE
      DOUBLE PRECISION SSMQCD
      COMPLEX ZI,ZONE,Z1(2),Z2(2)
      INTEGER IW,JZ,IZ,ISZJZ,ISWIW
      DATA FUDGE/1.0/,AMPI/.140/
      DATA ZONE/(1.,0.)/,ZI/(0.,1.)/,FPI/.1315/
C
C          Constants from neutralino mass matrix
C
      AMPL=2.4E18
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      E=SQRT(4*PI/128.)
C
      TANW=SQRT(SN2THW/(1.-SN2THW))
      COTW=1./TANW
      STHW=SQRT(SN2THW)
      CTHW=SQRT(1.-SN2THW)
      APL=.25*(3*TANW-COTW)
      BTL=.25*(COTW+TANW)
      APN=.25*(TANW+COTW)
      BTN=-.25*(COTW+TANW)
      APU=-5*TANW/12.+COTW/4.
      BTU=-.25*(COTW+TANW)
      APD=-COTW/4.+TANW/12.
      BTD=.25*(COTW+TANW)
C
      TANB=1./RV2V1
      BE=ATAN(1./RV2V1)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(AMTP))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(AMTP))
      FB=G*MBQ/SR2/AMW/COS(BE)
      FT=G*MTQ/SR2/AMW/SIN(BE)
      MTAMTA=AMTAU*(1.-SUALFE(AMTAU**2)/PI)
      MTAMB=MTAMTA*(SUALFE(AMBT**2)/SUALFE(AMTAU**2))**(-27./76.)
      MTAMZ=MTAMB*(SUALFE(AMZ**2)/SUALFE(AMBT**2))**(-27./80.)
      FL=G*MTAMZ/SR2/AMW/COS(BE)
      SNW1=SIGN(1.,AMW1SS)
      SNW2=SIGN(1.,AMW2SS)
      AMWISS(1)=AMW1SS
      AMWISS(2)=AMW2SS
      BWI(1)=-FT*SNW1*COS(GAMMAR)
      BWI(2)=FT*SNW2*THY*SIN(GAMMAR)
      BPWI(1)=-FB*COS(GAMMAL)
      BPWI(2)=FB*THX*SIN(GAMMAL)
      BPLWI(1)=-FL*COS(GAMMAL)
      BPLWI(2)=FL*THX*SIN(GAMMAL)
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
      MZ1=ABS(AMZ1SS)
      XWINO=.5*(THX*SIN(GAMMAL)*COS(GAMMAL)
     $-THY*SIN(GAMMAR)*COS(GAMMAR))
      YWINO=.5*(THX*SIN(GAMMAL)*COS(GAMMAL)
     $+THY*SIN(GAMMAR)*COS(GAMMAR))
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
C
C          Constants from Higgs mass matrix
C
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
C
C          Gaugino couplings
C
      DO 100 IZ=1,4
        AUI(IZ)=G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ)
        BUI(IZ)=4.*GP/3./SR2*ZMIXSS(4,IZ)
        ADI(IZ)=-G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ)
        BDI(IZ)=-2.*GP/3./SR2*ZMIXSS(4,IZ)
        ALI(IZ)=G/SR2*ZMIXSS(3,IZ)+GP/SR2*ZMIXSS(4,IZ)
        BLI(IZ)=-SR2*GP*ZMIXSS(4,IZ)
        ANI(IZ)=G/SR2*ZMIXSS(3,IZ)-GP/SR2*ZMIXSS(4,IZ)
        BNI(IZ)=0.0
100   CONTINUE
C
      DO 110 IZ=1,4
      DO 110 JZ=1,4
        IF(IZ.LT.JZ) THEN
          WIJ(IZ,JZ)=SQRT(G**2+GP**2)
     $    *(ZMIXSS(1,IZ)*ZMIXSS(1,JZ)-ZMIXSS(2,IZ)*ZMIXSS(2,JZ))/4.
        ELSEIF(IZ.GT.JZ) THEN
          WIJ(IZ,JZ)=-SQRT(G**2+GP**2)
     $    *(ZMIXSS(1,IZ)*ZMIXSS(1,JZ)-ZMIXSS(2,IZ)*ZMIXSS(2,JZ))/4.
        ENDIF
110   CONTINUE
C
      AUWI(2)=G*THX*COS(GAMMAL)
      ADWI(2)=SNW2*G*THY*COS(GAMMAR)
      ALWI(2)=ADWI(2)
      ANWI(2)=AUWI(2)
      AUWI(1)=G*SIN(GAMMAL)
      ADWI(1)=SNW1*G*SIN(GAMMAR)
      ALWI(1)=ADWI(1)
      ANWI(1)=AUWI(1)
C
      DO 120 IZ=1,4
        XIM(IZ)=.5*(SNW1*SIGN(1.,AMZISS(IZ))*(COS(GAMMAR)
     $  *ZMIXSS(1,IZ)/SR2+SIN(GAMMAR)*ZMIXSS(3,IZ))-COS(GAMMAL)
     $  *ZMIXSS(2,IZ)/SR2+SIN(GAMMAL)*ZMIXSS(3,IZ))
        YIM(IZ)=.5*(-SNW1*SIGN(1.,AMZISS(IZ))*(COS(GAMMAR)
     $  *ZMIXSS(1,IZ)/SR2+SIN(GAMMAR)*ZMIXSS(3,IZ))-COS(GAMMAL)
     $  *ZMIXSS(2,IZ)/SR2+SIN(GAMMAL)*ZMIXSS(3,IZ))
        XIP(IZ)=.5*(SNW2*SIGN(1.,AMZISS(IZ))*THY*(-SIN(GAMMAR)
     $  *ZMIXSS(1,IZ)/SR2+COS(GAMMAR)*ZMIXSS(3,IZ))+THX*(SIN(GAMMAL)
     $  *ZMIXSS(2,IZ)/SR2+COS(GAMMAL)*ZMIXSS(3,IZ)))
        YIP(IZ)=.5*(-SNW2*SIGN(1.,AMZISS(IZ))*THY*(-SIN(GAMMAR)
     $  *ZMIXSS(1,IZ)/SR2+COS(GAMMAR)*ZMIXSS(3,IZ))+THX*(SIN(GAMMAL)
     $  *ZMIXSS(2,IZ)/SR2+COS(GAMMAL)*ZMIXSS(3,IZ)))
120   CONTINUE
C
      DO 130 IZ=1,4
      DO 130 JZ=1,4
        IF(IZ.NE.JZ) THEN
          SNIJ(IZ,JZ)=-1.*SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
          XLIJ(IZ,JZ)=-SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
     $    *(ZMIXSS(2,IZ)*SINA-ZMIXSS(1,IZ)*COSA)
     $    *(G*ZMIXSS(3,JZ)-GP*ZMIXSS(4,JZ))/2.
          XHIJ(IZ,JZ)=-SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
     $    *(ZMIXSS(2,IZ)*COSA+ZMIXSS(1,IZ)*SINA)
     $    *(G*ZMIXSS(3,JZ)-GP*ZMIXSS(4,JZ))/2.
          XPIJ(IZ,JZ)=SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
     $    *(ZMIXSS(2,IZ)*SINBE-ZMIXSS(1,IZ)*COSBE)
     $    *(G*ZMIXSS(3,JZ)-GP*ZMIXSS(4,JZ))/2.
          HIJ(IZ,JZ)=-SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
     $    *(ZMIXSS(2,IZ)*COSA+ZMIXSS(1,IZ)*SINA)
     $    *(G*ZMIXSS(3,JZ)-GP*ZMIXSS(4,JZ))/2.
        ENDIF
130   CONTINUE
C
      SP=-.5*(-THY*SNW2*COSBE*SIN(GAMMAL)*SIN(GAMMAR)+
     $THY*SNW2*SINBE*COS(GAMMAL)*COS(GAMMAR)-
     $THX*SNW1*COSBE*COS(GAMMAL)*COS(GAMMAR)+
     $THX*SNW1*SINBE*SIN(GAMMAL)*SIN(GAMMAR))
      PP=-.5*(-THY*SNW2*COSBE*SIN(GAMMAL)*SIN(GAMMAR)+
     $THY*SNW2*SINBE*COS(GAMMAL)*COS(GAMMAR)+
     $THX*SNW1*COSBE*COS(GAMMAL)*COS(GAMMAR)-
     $THX*SNW1*SINBE*SIN(GAMMAL)*SIN(GAMMAR))
C
      SL=.5*(THY*SNW2*SINA*COS(GAMMAL)*COS(GAMMAR)-
     $THY*SNW2*COSA*SIN(GAMMAL)*SIN(GAMMAR)+
     $THX*SNW1*COSA*COS(GAMMAL)*COS(GAMMAR)-
     $THX*SNW1*SINA*SIN(GAMMAL)*SIN(GAMMAR))
      PL=.5*(THY*SNW2*SINA*COS(GAMMAL)*COS(GAMMAR)-
     $THY*SNW2*COSA*SIN(GAMMAL)*SIN(GAMMAR)-
     $THX*SNW1*COSA*COS(GAMMAL)*COS(GAMMAR)+
     $THX*SNW1*SINA*SIN(GAMMAL)*SIN(GAMMAR))
C
      SH=.5*(THY*SNW2*COSA*COS(GAMMAL)*COS(GAMMAR)+
     $THY*SNW2*SINA*SIN(GAMMAL)*SIN(GAMMAR)-
     $THX*SNW1*SINA*COS(GAMMAL)*COS(GAMMAR)-
     $THX*SNW1*COSA*SIN(GAMMAL)*SIN(GAMMAR))
      PH=.5*(THY*SNW2*COSA*COS(GAMMAL)*COS(GAMMAR)+
     $THY*SNW2*SINA*SIN(GAMMAL)*SIN(GAMMAR)+
     $THX*SNW1*SINA*COS(GAMMAL)*COS(GAMMAR)+
     $THX*SNW1*COSA*SIN(GAMMAL)*SIN(GAMMAR))
C
      DO 140 IZ=1,4
        V1I(IZ)=-SIN(GAMMAR)/SR2*(G*ZMIXSS(3,IZ)+GP*ZMIXSS(4,IZ))
     $  -G*COS(GAMMAR)*ZMIXSS(1,IZ)
        V2I(IZ)=COS(GAMMAR)/SR2*(G*ZMIXSS(3,IZ)+GP*ZMIXSS(4,IZ))
     $  -G*SIN(GAMMAR)*ZMIXSS(1,IZ)
        V3I(IZ)=-SIN(GAMMAL)/SR2*(G*ZMIXSS(3,IZ)+GP*ZMIXSS(4,IZ))
     $  +G*COS(GAMMAL)*ZMIXSS(2,IZ)
        V4I(IZ)=COS(GAMMAL)/SR2*(G*ZMIXSS(3,IZ)+GP*ZMIXSS(4,IZ))
     $  +G*SIN(GAMMAL)*ZMIXSS(2,IZ)
140   CONTINUE
C
C-----------------------------------------------------------------------
C          Generate Chargino Branching Fractions
C-----------------------------------------------------------------------
      ISZ(1)=ISZ1
      ISZ(2)=ISZ2
      ISZ(3)=ISZ3
      ISZ(4)=ISZ4
C     FIRST TRY EXCLUSIVE DECAY TO SINGLE PION
      IF (MW1.GT.(MZ1+AMPI).AND.MW1.LT.(MZ1+1.)) THEN
        WID=G**4*FPI**2*SQRT(SSXLAM(MW1**2,MZ1**2,AMPI**2))/
     $      128./MW1**3/PI/AMW**4*((XIM(1)**2+YIM(1)**2)*
     $      (MW1-MZ1)**2*(MW1+MZ1)**2-AMPI**2*(XIM(1)**2*(MW1-MZ1)**2+
     $      YIM(1)**2*(MW1+MZ1)**2))
        CALL SSSAVE(ISW1,WID,ISZ(1),120,0,0,0)
      END IF
      DO 300 IW=1,2
C          Loop over w1, w2
        IF(IW.EQ.1) THEN
            MWIW=MW1
            SNIW=SNW1
            ISWIW=ISW1
        ELSE
            MWIW=MW2
            SNIW=SNW2
            ISWIW=ISW2
        ENDIF
C
C          Decays to zj
C
        DO 310 JZ=1,4
          MZJZ=ABS(AMZISS(JZ))
          ISZJZ=ISZ(JZ)
          SNJZ=SIGN(1.,AMZISS(JZ))
          THJZ=0.
          IF (AMZISS(JZ).LT.0.) THJZ=1.
          IF(MWIW.LE.FUDGE*MZJZ) GOTO 310
C          Couplings
          IF(IW.EQ.1) THEN
            XIPM=XIM(JZ)
            YIPM=YIM(JZ)
          ELSE
            XIPM=XIP(JZ)
            YIPM=YIP(JZ)
          ENDIF
C
C          wi --> f + fbar + zj
C
          IF (MWIW.GT.(MZJZ+AMUP+AMDN)) THEN
          IF (IW.EQ.1.AND.JZ.EQ.1.AND.(MW1.LT.(MZ1+1.))) THEN
C           EXIT HADRONIC MODE IF EXCLUSIVE PION DECAY ACTIVE
            GO TO 200
          END IF
            IF (MWIW.LT.(AMW+MZJZ)) THEN
            TMP(1)=XIPM**2+YIPM**2
            TMP(2)=XIPM**2-YIPM**2
            TMP(3)=MWIW
            TMP(4)=MZJZ
            TERM1=SSXINT(MZJZ,SSWZF1,(MWIW**2+MZJZ**2)/2./MWIW)/
     $            2./MWIW/(2*PI)**5*2*G**4*PI**2/3.
            ELSE
            TERM1=0.
            END IF
            ULIM=MWIW/2.*(1.-MZJZ**2/MWIW**2)
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.AMULSS) THEN
            TMP(2)=AMULSS
            PSIINT=SSXINT(0.,SSWZF2,ULIM)
            T1=AUI(JZ)**2*ADWI(IW)**2*PSIINT
            ELSE
            T1=0.
            END IF
            IF (MWIW.LT.AMDLSS) THEN
            TMP(2)=AMDLSS
            PSIINT=SSXINT(0.,SSWZF2,ULIM)
            T2=ADI(JZ)**2*AUWI(IW)**2*PSIINT
            ELSE
            T2=0.
            END IF
            IF (MWIW.LT.AMDLSS.AND.MWIW.LT.AMULSS) THEN
            TMP(2)=0.
            TMP(4)=AMDLSS
            TMP(5)=AMULSS
            PHIINT=SSXINT(0.,SSGX2,ULIM)
            CC=2*SIGN(1.,AMZISS(JZ))*AUWI(IW)*ADWI(IW)*AUI(JZ)*ADI(JZ)
            T3=CC*PHIINT
            ELSE
            T3=0.
            END IF
            TERM2=(T1+T2+T3)/2./MWIW/(2*PI)**5
            FACTOR=1./2./MWIW/(2*PI)**5*2*SR2*G**2
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMULSS) THEN
            TMP(2)=AMULSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM3=FACTOR*ADWI(IW)*AUI(JZ)*((XIPM-YIPM)*XI1
     $      -(XIPM+YIPM)*XI2)*SIGN(1.,AMZISS(JZ))
            ELSE
            TERM3=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMDLSS) THEN
            TMP(2)=AMDLSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM4=-FACTOR*AUWI(IW)*ADI(JZ)*((XIPM+YIPM)*XI1
     $      -(XIPM-YIPM)*XI2)
            ELSE
            TERM4=0.
            END IF
            WID=3*(TERM1+TERM2+TERM3+TERM4)
            CALL SSSAVE(ISWIW,WID,ISZJZ,IDUP,-IDDN,0,0)
C         Enter information for decay matrix element
            Z1(1)=ZI**THJZ*G*XIPM
            Z1(2)=ZI**THJZ*G*YIPM
            Z2(1)=G/2./SR2
            Z2(2)=-G/2./SR2
            CALL SSME3(1,AMW,Z1,Z2)
            Z1(1)=ZI*AUWI(IW)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ))/2.
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMDLSS,Z1,Z2)
            Z1(1)=CONJG(ZI*ADWI(IW))/2.
            Z1(2)=-Z1(1)
            Z2(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMULSS,Z1,Z2)
          END IF
C          wi --> c + sbar + zj
200       IF (MWIW.GT.(MZJZ+AMCH+AMST)) THEN
            IF (MWIW.LT.(AMW+MZJZ)) THEN
            TMP(1)=XIPM**2+YIPM**2
            TMP(2)=XIPM**2-YIPM**2
            TMP(3)=MWIW
            TMP(4)=MZJZ
            TERM1=SSXINT(MZJZ,SSWZF1,(MWIW**2+MZJZ**2)/2./MWIW)/
     $            2./MWIW/(2*PI)**5*2*G**4*PI**2/3.
            ELSE
            TERM1=0.
            END IF
            ULIM=MWIW/2.*(1.-MZJZ**2/MWIW**2)
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.AMCLSS) THEN
            TMP(2)=AMCLSS
            PSIINT=SSXINT(0.,SSWZF2,ULIM)
            T1=AUI(JZ)**2*ADWI(IW)**2*PSIINT
            ELSE
            T1=0.
            END IF
            IF (MWIW.LT.AMSLSS) THEN
            TMP(2)=AMSLSS
            PSIINT=SSXINT(0.,SSWZF2,ULIM)
            T2=ADI(JZ)**2*AUWI(IW)**2*PSIINT
            ELSE
            T2=0.
            END IF
            IF (MWIW.LT.AMSLSS.AND.MWIW.LT.AMCLSS) THEN
            TMP(2)=0.
            TMP(4)=AMSLSS
            TMP(5)=AMCLSS
            PHIINT=SSXINT(0.,SSGX2,ULIM)
            CC=2*SIGN(1.,AMZISS(JZ))*AUWI(IW)*ADWI(IW)*AUI(JZ)*ADI(JZ)
            T3=CC*PHIINT
            ELSE
            T3=0.
            END IF
            TERM2=(T1+T2+T3)/2./MWIW/(2*PI)**5
            FACTOR=1./2./MWIW/(2*PI)**5*2*SR2*G**2
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMCLSS) THEN
            TMP(2)=AMCLSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM3=FACTOR*ADWI(IW)*AUI(JZ)*((XIPM-YIPM)*XI1
     $      -(XIPM+YIPM)*XI2)*SIGN(1.,AMZISS(JZ))
            ELSE
            TERM3=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMSLSS) THEN
            TMP(2)=AMSLSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM4=-FACTOR*AUWI(IW)*ADI(JZ)*((XIPM+YIPM)*XI1
     $      -(XIPM-YIPM)*XI2)
            ELSE
            TERM4=0.
            END IF
            WID=3*(TERM1+TERM2+TERM3+TERM4)
            CALL SSSAVE(ISWIW,WID,ISZJZ,IDCH,-IDST,0,0)
C         Enter information for decay matrix element
            Z1(1)=ZI**THJZ*G*XIPM
            Z1(2)=ZI**THJZ*G*YIPM
            Z2(1)=G/2./SR2
            Z2(2)=-G/2./SR2
            CALL SSME3(1,AMW,Z1,Z2)
            Z1(1)=ZI*AUWI(IW)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ))/2.
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMSLSS,Z1,Z2)
            Z1(1)=CONJG(ZI*ADWI(IW))/2.
            Z1(2)=-Z1(1)
            Z2(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMCLSS,Z1,Z2)
          END IF
C          wi -> t + bbar + zj neglected since 2-body modes should dominate
C          wi --> nu_e + e + zj
          IF (MWIW.GT.(MZJZ+AME)) THEN
            IF (MWIW.LT.(AMW+MZJZ)) THEN
            TMP(1)=XIPM**2+YIPM**2
            TMP(2)=XIPM**2-YIPM**2
            TMP(3)=MWIW
            TMP(4)=MZJZ
            TERM1=SSXINT(MZJZ,SSWZF1,(MWIW**2+MZJZ**2)/2./MWIW)/
     $            2./MWIW/(2*PI)**5*2*G**4*PI**2/3.
            ELSE
            TERM1=0.
            END IF
            ULIM=MWIW/2.*(1.-MZJZ**2/MWIW**2)
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.AMN1SS) THEN
            TMP(2)=AMN1SS
            T1=ANI(JZ)**2*ALWI(IW)**2*SSXINT(0.,SSWZF2,ULIM)
            ELSE
            T1=0.
            END IF
            IF (MWIW.LT.AMELSS) THEN
            TMP(2)=AMELSS
            T2=ALI(JZ)**2*ANWI(IW)**2*SSXINT(0.,SSWZF2,ULIM)
            ELSE
            T2=0.
            END IF
            IF (MWIW.LT.AMELSS.AND.MWIW.LT.AMN1SS) THEN
            TMP(2)=0.
            TMP(4)=AMELSS
            TMP(5)=AMN1SS
            T3=-2*SIGN(1.,AMZISS(JZ))*ANWI(IW)*ALWI(IW)*ANI(JZ)*
     $         ALI(JZ)*SSXINT(0.,SSGX2,ULIM)
            ELSE
            T3=0.
            END IF
            TERM2=(T1+T2+T3)/2./MWIW/(2*PI)**5
            FACTOR=1./2./MWIW/(2*PI)**5*2*SR2*G**2
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMN1SS) THEN
            TMP(2)=AMN1SS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM3=FACTOR*ALWI(IW)*ANI(JZ)*((XIPM-YIPM)*XI1
     $      -(XIPM+YIPM)*XI2)*SIGN(1.,AMZISS(JZ))
            ELSE
            TERM3=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMELSS) THEN
            TMP(2)=AMELSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM4=FACTOR*ANWI(IW)*ALI(JZ)*((XIPM+YIPM)*XI1
     $      -(XIPM-YIPM)*XI2)
            ELSE
            TERM4=0.
            END IF
            WID=TERM1+TERM2+TERM3+TERM4
            CALL SSSAVE(ISWIW,WID,ISZJZ,-IDE,IDNE,0,0)
C         Enter information for decay matrix element
            Z1(1)=ZI**THJZ*G*XIPM
            Z1(2)=ZI**THJZ*G*YIPM
            Z2(1)=G/2./SR2
            Z2(2)=-G/2./SR2
            CALL SSME3(1,AMW,Z1,Z2)
            Z1(1)=ZI*ANWI(IW)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ))/2.
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMELSS,Z1,Z2)
            Z1(1)=CONJG(ZI*ALWI(IW))/2.
            Z1(2)=-Z1(1)
            Z2(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMN1SS,Z1,Z2)
          END IF
C          wi --> nu_mu + mu + zj
          IF (MWIW.GT.(MZJZ+AMMU)) THEN
            IF (MWIW.LT.(AMW+MZJZ)) THEN
            TMP(1)=XIPM**2+YIPM**2
            TMP(2)=XIPM**2-YIPM**2
            TMP(3)=MWIW
            TMP(4)=MZJZ
            TERM1=SSXINT(MZJZ,SSWZF1,(MWIW**2+MZJZ**2)/2./MWIW)/
     $            2./MWIW/(2*PI)**5*2*G**4*PI**2/3.
            ELSE
            TERM1=0.
            END IF
            ULIM=MWIW/2.*(1.-MZJZ**2/MWIW**2)
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.AMN2SS) THEN
            TMP(2)=AMN2SS
            T1=ANI(JZ)**2*ALWI(IW)**2*SSXINT(0.,SSWZF2,ULIM)
            ELSE
            T1=0.
            END IF
            IF (MWIW.LT.AMMLSS) THEN
            TMP(2)=AMMLSS
            T2=ALI(JZ)**2*ANWI(IW)**2*SSXINT(0.,SSWZF2,ULIM)
            ELSE
            T2=0.
            END IF
            IF (MWIW.LT.AMMLSS.AND.MWIW.LT.AMN2SS) THEN
            TMP(2)=0.
            TMP(4)=AMMLSS
            TMP(5)=AMN2SS
            T3=-2*SIGN(1.,AMZISS(JZ))*ANWI(IW)*ALWI(IW)*ANI(JZ)*
     $         ALI(JZ)*SSXINT(0.,SSGX2,ULIM)
            ELSE
            T3=0.
            END IF
            TERM2=(T1+T2+T3)/2./MWIW/(2*PI)**5
            FACTOR=1./2./MWIW/(2*PI)**5*2*SR2*G**2
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMN2SS) THEN
            TMP(2)=AMN2SS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM3=FACTOR*ALWI(IW)*ANI(JZ)*((XIPM-YIPM)*XI1
     $      -(XIPM+YIPM)*XI2)*SIGN(1.,AMZISS(JZ))
            ELSE
            TERM3=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMMLSS) THEN
            TMP(2)=AMMLSS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERM4=FACTOR*ANWI(IW)*ALI(JZ)*((XIPM+YIPM)*XI1
     $      -(XIPM-YIPM)*XI2)
            ELSE
            TERM4=0.
            END IF
            WID=TERM1+TERM2+TERM3+TERM4
            CALL SSSAVE(ISWIW,WID,ISZJZ,-IDMU,IDNM,0,0)
C         Enter information for decay matrix element
            Z1(1)=ZI**THJZ*G*XIPM
            Z1(2)=ZI**THJZ*G*YIPM
            Z2(1)=G/2./SR2
            Z2(2)=-G/2./SR2
            CALL SSME3(1,AMW,Z1,Z2)
            Z1(1)=ZI*ANWI(IW)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ))/2.
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMMLSS,Z1,Z2)
            Z1(1)=CONJG(ZI*ALWI(IW))/2.
            Z1(2)=-Z1(1)
            Z2(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.
            Z2(2)=Z2(1)
            CALL SSME3(3,AMN2SS,Z1,Z2)
          END IF
C          wi --> nu_tau + tau + zj ; includes mixing and Yukawas
          FACT=1./2./MWIW/(2*PI)**5
          ALJZ1=-ALI(JZ)*COSL-FL*ZMIXSS(2,JZ)*SINL
          ALJZ2=-ALI(JZ)*SINL+FL*ZMIXSS(2,JZ)*COSL
          BEJZ1=BLI(JZ)*SINL+FL*ZMIXSS(2,JZ)*COSL
          BEJZ2=-BLI(JZ)*COSL+FL*ZMIXSS(2,JZ)*SINL
          SNJZ=SIGN(1.,AMZISS(JZ))
C         Change ALWI def'ns in accord with Drees note
          ALWI(1)=-G*SIN(GAMMAR)
          ALWI(2)=-G*THY*COS(GAMMAR)
C          Polarization for stau_1 -> z1ss+tau.
C          See above for other cases.
          IF(IW.EQ.1.AND.JZ.EQ.1) THEN
            PTAU1(JZ)=(BEJZ1**2-ALJZ1**2)/(BEJZ1**2+ALJZ1**2)
            PTAU2(JZ)=(BEJZ2**2-ALJZ2**2)/(BEJZ2**2+ALJZ2**2)
          ENDIF
          IF (IW.EQ.1) THEN
            ALIW1=-G*SIN(GAMMAL)*COSL+FL*COS(GAMMAL)*SINL
            ALIW2=-G*SIN(GAMMAL)*SINL-FL*COS(GAMMAL)*COSL
            AHCJZ=COSBE*V2I(JZ)
            BHCJZ=-SINBE*V4I(JZ)
          ELSE IF (IW.EQ.2) THEN
            ALIW1=(-G*COS(GAMMAL)*COSL-FL*SIN(GAMMAL)*SINL)*THX
            ALIW2=(-G*COS(GAMMAL)*SINL+FL*SIN(GAMMAL)*COSL)*THX
            AHCJZ=COSBE*V1I(JZ)*THY
            BHCJZ=-SINBE*V3I(JZ)*THX
          END IF
          IF (MWIW.GT.(MZJZ+AMTAU)) THEN
            IF (MWIW.LT.(AMW+MZJZ)) THEN
            TMP(1)=XIPM**2+YIPM**2
            TMP(2)=XIPM**2-YIPM**2
            TMP(3)=MWIW
            TMP(4)=MZJZ
            TERMW=SSXINT(MZJZ,SSWZF1,(MWIW**2+MZJZ**2)/2./MWIW)*
     $            FACT*2*G**4*PI**2/3.
            ELSE
            TERMW=0.
            END IF
            ULIM=(MWIW**2-MZJZ**2)/2./MWIW
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.AMN3SS) THEN
            TMP(2)=AMN3SS
            POLNL=FACT*ANI(JZ)**2*ALWI(IW)**2*
     $            SSXINT(0.,SSWZF2,ULIM)
            POLNR=POLNL*BPLWI(IW)**2/ALWI(IW)**2
            TERMN=POLNL+POLNR
            ELSE
            POLNL=0.
            POLNR=0.
            TERMN=0.
            END IF
            IF (MWIW.LT.AML1SS) THEN
            TMP(2)=AML1SS
            POL1L=FACT*ALJZ1**2*ALIW1**2*SSXINT(0.,SSWZF2,ULIM)
            POL1R=POL1L*BEJZ1**2/ALJZ1**2
            TERM1=POL1L+POL1R
            ELSE
            POL1L=0.
            POL1R=0.
            TERM1=0.
            END IF
            IF (MWIW.LT.AML2SS) THEN
            TMP(2)=AML2SS
            POL2L=FACT*ALJZ2**2*ALIW2**2*SSXINT(0.,SSWZF2,ULIM)
            POL2R=POL2L*BEJZ2**2/ALJZ2**2
            TERM2=POL2L+POL2R
            ELSE
            POL2L=0.
            POL2R=0.
            TERM2=0.
            END IF
            IF (MWIW.LT.AML1SS) THEN
            TMP(2)=0.
            TMP(3)=MZJZ
            TMP(4)=AML1SS
            TMP(5)=AML2SS
            POL12L=FACT*2*ALIW1*ALIW2*ALJZ1*ALJZ2*
     $             SSXINT(0.,SSGX1,ULIM)
            POL12R=POL12L*BEJZ1*BEJZ2/ALJZ1/ALJZ2
            TERM12=POL12L+POL12R
            ELSE
            POL12L=0.
            POL12R=0.
            TERM12=0.
            END IF
            IF (MWIW.LT.(AMHC+MZJZ)) THEN
            TMP(2)=AMHC
            TMP(3)=MZJZ
            TMP(4)=AHCJZ
            TMP(5)=BHCJZ
            TMP(6)=SIGN(1.,AMZISS(JZ))*SIGN(1.,AMWISS(IW))
            TERMH=FACT*PI**2*MWIW*(G*MTAMZ*TANB/AMW)**2/2.*
     $         SSXINT(MZJZ,SSWZF6,(MWIW**2+MZJZ**2)/2./MWIW)
            ELSE
            TERMH=0.
            END IF
            IF (MWIW.LT.AML1SS.AND.MWIW.LT.AMN3SS) THEN
            TMP(2)=0.
            TMP(3)=MZJZ
            TMP(4)=AMN3SS
            TMP(5)=AML1SS
            POLN1L=+2*FACT*ANI(JZ)*ALIW1*SNJZ*SNIW*ALWI(IW)*
     $             ALJZ1*SSXINT(0.,SSGX2,ULIM)
            POLN1R=-2*FACT*ANI(JZ)*ALIW1*BPLWI(IW)*BEJZ1*
     $             SSXINT(0.,SSGX8,ULIM)
            TERMN1=POLN1L+POLN1R
            ELSE
            POLN1L=0.
            POLN1R=0.
            TERMN1=0.
            END IF
            IF (MWIW.LT.AML2SS.AND.MWIW.LT.AMN3SS) THEN
            TMP(2)=0.
            TMP(3)=MZJZ
            TMP(4)=AMN3SS
            TMP(5)=AML2SS
            POLN2L=+2*FACT*ANI(JZ)*ALIW2*SNJZ*SNIW*ALWI(IW)*ALJZ2*
     $             SSXINT(0.,SSGX2,ULIM)
            POLN2R=-2*FACT*ANI(JZ)*ALIW2*BPLWI(IW)*BEJZ2*
     $             SSXINT(0.,SSGX8,ULIM)
            TERMN2=POLN2L+POLN2R
            ELSE
            POLN2L=0.
            POLN2R=0.
            TERMN2=0.
            END IF
            TMP(1)=MWIW
            TMP(3)=MZJZ
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AMN3SS) THEN
            TMP(2)=AMN3SS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERMWN=2*SR2*G**2*FACT*ALWI(IW)*ANI(JZ)*((XIPM-
     $         YIPM)*XI1-(XIPM+YIPM)*XI2)*SIGN(1.,AMZISS(JZ))
            ELSE
            TERMWN=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AML1SS) THEN
            TMP(2)=AML1SS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERMW1=2*SR2*G**2*FACT*ALIW1*ALJZ1*
     $             ((XIPM+YIPM)*XI1-(XIPM-YIPM)*XI2)
            ELSE
            TERMW1=0.
            END IF
            IF (MWIW.LT.(MZJZ+AMW).AND.MWIW.LT.AML2SS) THEN
            TMP(2)=AML2SS
            XI1=SSXINT(0.,SSWZF4,(MWIW-MZJZ)**2)
            XI2=SSXINT(0.,SSWZF5,(MWIW-MZJZ)**2)
            TERMW2=2*SR2*G**2*FACT*ALIW2*ALJZ2*
     $             ((XIPM+YIPM)*XI1-(XIPM-YIPM)*XI2)
            ELSE
            TERMW2=0.
            END IF
            TMP(2)=MZJZ
            TMP(3)=AMHC
            TMP(5)=AHCJZ
            TMP(6)=BHCJZ
            TMP(7)=SNJZ*SNIW
            IF (MWIW.LT.(AMHC+MZJZ).AND.MWIW.LT.AML1SS) THEN
            TMP(4)=AML1SS
            TERMH1=PI**2/2./MWIW*FACT*SR2*ALIW1*BEJZ1*G*MTAMZ*
     $             TANB/AMW*SSXINT(0.,SSWZF7,(MWIW-MZJZ)**2)
            ELSE
            TERMH1=0.
            END IF
            IF (MWIW.LT.(AMHC+MZJZ).AND.MWIW.LT.AML2SS) THEN
            TMP(4)=AML2SS
            TERMH2=PI**2/2./MWIW*FACT*SR2*ALIW2*BEJZ2*G*MTAMZ*
     $             TANB/AMW*SSXINT(0.,SSWZF7,(MWIW-MZJZ)**2)
            ELSE
            TERMH2=0.
            END IF
            IF (MWIW.LT.(AMHC+MZJZ).AND.MWIW.LT.AMN3SS) THEN
            TMP(4)=AMN3SS
            TERMHN=PI**2/2./MWIW*FACT*SR2*ANI(JZ)*BPLWI(IW)*G*
     $         MTAMZ*TANB/AMW*SSXINT(0.,SSWZF7,(MWIW-MZJZ)**2)
            ELSE
            TERMHN=0.
            END IF
            WID=TERMW+TERMN+TERM1+TERM2+TERMH+TERMWN+TERMW1+
     $          TERMW2+TERM12+TERMN1+TERMN2+TERMH1+TERMH2+
     $          TERMHN
C          tau polarization for 3-body w1 -> z1 tau nu
            IF (IW.EQ.1.AND.JZ.EQ.1.AND.WID.GT.0.) THEN
              PTAUWZ=(POLNR+POL1R+POL2R+POL12R+TERMH+POLN1R+
     $               POLN2R+TERMHN+TERMH1+TERMH2-(TERMW+POLNL+
     $               POL1L+POL2L+POL12L+POLN1L+POLN2L+TERMWN+
     $               TERMW1+TERMW2))/WID
            END IF
            CALL SSSAVE(ISWIW,WID,ISZJZ,-IDTAU,IDNT,0,0)
            Z1(1)=ZI**THJZ*G*XIPM
            Z1(2)=ZI**THJZ*G*YIPM
            Z2(1)=-G/2./SR2
            Z2(2)=-Z2(1)
            CALL SSME3(1,AMW,Z1,Z2)
            Z1(1)=-(ZI)**THJZ*ALIW1/2.
            Z1(2)=-Z1(1)
            Z2(1)=(ALJZ1+(-1.)**THJZ*BEJZ1)/2.
            Z2(2)=(ALJZ1-(-1.)**THJZ*BEJZ1)/2.
            CALL SSME3(2,AML1SS,Z1,Z2)
            Z1(1)=-(ZI)**THJZ*ALIW2/2.
            Z1(2)=-Z1(1)
            Z2(1)=(ALJZ2+(-1.)**THJZ*BEJZ2)/2.
            Z2(2)=(ALJZ2-(-1.)**THJZ*BEJZ2)/2.
            CALL SSME3(2,AML2SS,Z1,Z2)
            Z1(1)=(-1.)**THJZ*ANI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=(SNIW*ALWI(IW)+BPLWI(IW))/2.
            Z2(2)=(SNIW*ALWI(IW)-BPLWI(IW))/2.
            CALL SSME3(3,AMN3SS,Z1,Z2)
            Z1(1)=ZI**THJZ*G*MTAMZ*TANB/SR2/AMW/2.
            Z1(2)=-Z1(1)
            Z2(1)=(SNIW*AHCJZ-SNJZ*BHCJZ)/2.
            Z2(2)=(SNIW*AHCJZ+SNJZ*BHCJZ)/2.
            CALL SSME3(4,AMHC,Z1,Z2)
          END IF
C
C          wi --> w + zj
C
          IF (MWIW.GT.(MZJZ+AMW)) THEN
            EF=MWIW**2+MZJZ**2-AMW**2+((MWIW**2-MZJZ**2)**2
     $         -AMW**4)/AMW/AMW
            WID=G*G*SQRT(SSXLAM(MWIW**2,MZJZ**2,AMW**2))/32./PI/
     $          MWIW**3*(2.*EF*(XIPM**2+YIPM**2)-12*MWIW*MZJZ*
     $          (XIPM**2-YIPM**2))
            CALL SSSAVE(ISWIW,WID,ISZJZ,IDW,0,0,0)
          END IF
C
C          wi --> h+ + zj
C
          IF (MWIW.GT.(MZJZ+AMHC)) THEN
            IF(IW.EQ.1) THEN
              A=(SNW1*COSBE*V2I(JZ)
     $        -SIGN(1.,AMZISS(JZ))*SINBE*V4I(JZ))/2.
              B=(SNW1*COSBE*V2I(JZ)
     $        +SIGN(1.,AMZISS(JZ))*SINBE*V4I(JZ))/2.
            ELSE
              A=(THY*SNW2*COSBE*V1I(JZ)
     $        -THX*SIGN(1.,AMZISS(JZ))*SINBE*V3I(JZ))/2.
              B=(THY*SNW2*COSBE*V1I(JZ)
     $        +THX*SIGN(1.,AMZISS(JZ))*SINBE*V3I(JZ))/2.
            ENDIF
            WID=SQRT(MWIW**4+MZJZ**4+AMHC**4-2.*(MWIW*MZJZ)**2-
     $      2.*(MWIW*AMHC)**2-2.*(MZJZ*AMHC)**2)/8./PI/MWIW**3*
     $      ((A*A+B*B)*(MWIW*MWIW+MZJZ*MZJZ-AMHC*AMHC)/2.
     $      +(A*A-B*B)*MWIW*MZJZ) 
            CALL SSSAVE(ISWIW,WID,ISZJZ,ISHC,0,0,0)
          ENDIF
310     CONTINUE
C
C          wi --> quark + squark
C
        IF(MWIW.GT.(AMULSS+AMDN)) THEN
          WID=3.*ADWI(IW)**2/32./PI/MWIW*(1.+AMDN**2/MWIW**2-
     $     AMULSS**2/MWIW**2)*SQRT(SSXLAM(MWIW**2,AMDN**2,AMULSS**2))
          CALL SSSAVE(ISWIW,WID,+ISUPL,-IDDN,0,0,0)
        END IF
        IF(MWIW.GT.(AMDLSS+AMUP)) THEN
          WID=3.*AUWI(IW)**2/32./PI/MWIW*(1.+AMUP**2/MWIW**2-
     $     AMDLSS**2/MWIW**2)*SQRT(SSXLAM(MWIW**2,AMUP**2,AMDLSS**2))
          CALL SSSAVE(ISWIW,WID,-ISDNL,+IDUP,0,0,0)
        END IF
        IF(MWIW.GT.(AMCLSS+AMST)) THEN
          WID=3.*ADWI(IW)**2/32./PI/MWIW*(1.+AMST**2/MWIW**2-
     $     AMCLSS**2/MWIW**2)*SQRT(SSXLAM(MWIW**2,AMST**2,AMCLSS**2))
          CALL SSSAVE(ISWIW,WID,+ISCHL,-IDST,0,0,0)
        END IF
        IF(MWIW.GT.(AMSLSS+AMCH)) THEN
          WID=3.*AUWI(IW)**2/32./PI/MWIW*(1.+AMCH**2/MWIW**2-
     $     AMCLSS**2/MWIW**2)*SQRT(SSXLAM(MWIW**2,AMCH**2,AMCLSS**2))
          CALL SSSAVE(ISWIW,WID,-ISSTL,+IDCH,0,0,0)
        ENDIF
        IF(MWIW.GT.(AMT1SS+AMBT)) THEN
          AS=(-ADWI(IW)*COST-BWI(IW)*SINT)**2
          WID=3*((AS+BPWI(IW)**2*COST**2)*(MWIW**2+AMBT**2-AMT1SS**2)
     $        +4*SQRT(AS)*BPWI(IW)*COST*AMBT*MWIW)/32./PI/MWIW**3*
     $      SQRT(SSXLAM(MWIW**2,AMBT**2,AMT1SS**2))
          CALL SSSAVE(ISWIW,WID,ISTP1,-IDBT,0,0,0)
        ENDIF
        IF(MWIW.GT.(AMT2SS+AMBT)) THEN
          AS=(-ADWI(IW)*SINT+BWI(IW)*COST)**2
          WID=3*((AS+BPWI(IW)**2*SINT**2)*(MWIW**2+AMBT**2-AMT2SS**2)
     $        +4*SQRT(AS)*BPWI(IW)*SINT*AMBT*MWIW)/32./PI/MWIW**3*
     $      SQRT(SSXLAM(MWIW**2,AMBT**2,AMT2SS**2))
          CALL SSSAVE(ISWIW,WID,ISTP2,-IDBT,0,0,0)
        ENDIF
        IF(MWIW.GT.(AMB1SS+AMTP)) THEN
          AS=(-AUWI(IW)*COSB-BPWI(IW)*SINB)**2
          WID=3*((AS+BWI(IW)**2*COSB**2)*(MWIW**2+AMTP**2-AMB1SS**2)
     $        +4*SQRT(AS)*BWI(IW)*COSB*AMTP*MWIW)/32./PI/MWIW**3*
     $      SQRT(SSXLAM(MWIW**2,AMTP**2,AMB1SS**2))
          CALL SSSAVE(ISWIW,WID,-ISBT1,IDTP,0,0,0)
        ENDIF
        IF(MWIW.GT.(AMB2SS+AMTP)) THEN
          AS=(-AUWI(IW)*SINB+BPWI(IW)*COSB)**2
          WID=3*((AS+BWI(IW)**2*SINB**2)*(MWIW**2+AMTP**2-AMB2SS**2)
     $        +4*SQRT(AS)*BWI(IW)*SINB*AMTP*MWIW)/32./PI/MWIW**3*
     $      SQRT(SSXLAM(MWIW**2,AMTP**2,AMB2SS**2))
          CALL SSSAVE(ISWIW,WID,-ISBT2,IDTP,0,0,0)
        ENDIF
C
C          wi --> lepton + slepton
C
        IF(MWIW.GT.(AMN1SS+AME)) THEN
          AS=(-ALWI(IW))**2
          WID=AS*(MWIW**2+AME**2-AMN1SS**2)/32./PI/MWIW**3*
     $        SQRT(SSXLAM(MWIW**2,AME**2,AMN1SS**2))
          CALL SSSAVE(ISWIW,WID,ISNEL,-IDE,0,0,0)
        END IF
        IF(MWIW.GT.(AMN2SS+AMMU)) THEN
          AS=(-ALWI(IW))**2
          WID=AS*(MWIW**2+AMMU**2-AMN2SS**2)/32./PI/MWIW**3*
     $        SQRT(SSXLAM(MWIW**2,AMMU**2,AMN2SS**2))
          CALL SSSAVE(ISWIW,WID,ISNML,-IDMU,0,0,0)
        END IF
        IF(MWIW.GT.(AMN3SS+AMTAU)) THEN
          AS=(-ALWI(IW))**2
          WID=((AS+BPLWI(IW)**2)*(MWIW**2+AMTAU**2-AMN3SS**2)+
     $        4*SQRT(AS)*BPLWI(IW)*AMTAU*MWIW)/32./PI/MWIW**3*
     $        SQRT(SSXLAM(MWIW**2,AMTAU**2,AMN3SS**2))
          CALL SSSAVE(ISWIW,WID,ISNTL,-IDTAU,0,0,0)
        ENDIF
        IF(MWIW.GT.AMELSS) THEN
          WID=ANWI(IW)**2*(MWIW**2-AMELSS**2)**2/32./PI/MWIW**3
          CALL SSSAVE(ISWIW,WID,-ISEL,IDNE,0,0,0)
        ENDIF
        IF(MWIW.GT.AMMLSS) THEN
          WID=ANWI(IW)**2*(MWIW**2-AMMLSS**2)**2/32./PI/MWIW**3
          CALL SSSAVE(ISWIW,WID,-ISMUL,IDNM,0,0,0)
        END IF
        IF(MWIW.GT.AML1SS) THEN
          AS=(-ANWI(IW)*COSL-BPLWI(IW)*SINL)**2
          WID=AS*(MWIW**2-AML1SS**2)**2/32./PI/MWIW**3
          CALL SSSAVE(ISWIW,WID,-ISTAU1,IDNT,0,0,0)
        END IF
        IF(MWIW.GT.AML2SS) THEN
          AS=(-ANWI(IW)*SINL+BPLWI(IW)*COSL)**2
          WID=AS*(MWIW**2-AML2SS**2)**2/32./PI/MWIW**3
          CALL SSSAVE(ISWIW,WID,-ISTAU2,IDNT,0,0,0)
        END IF
300   CONTINUE
C
C          w2 --> w1 + z
C          w2 --> w1 + f + fbar
C
      IF (MW2.GT.(MW1+AMZ)) THEN
        EF=MW2**2+MW1**2-AMZ**2+((MW2**2-MW1**2)**2-
     $     AMZ**4)/AMZ/AMZ
        Y=(THX*SIN(GAMMAL)*COS(GAMMAL)-THY*SIN(GAMMAR)*COS(GAMMAR))/2.
        Z=(THX*SIN(GAMMAL)*COS(GAMMAL)+THY*SIN(GAMMAR)*COS(GAMMAR))/2.
        WID=(COTW+TANW)**2*SQRT(SSXLAM(MW2**2,MW1**2,AMZ**2))/
     $      128./32./MW2**3*(2*EF*(Y*Y+Z*Z)+
     $      12*MW2*MW1*(Y*Y-Z*Z)*SNW2*SNW1)
        CALL SSSAVE(ISW2,WID,ISW1,IDZ,0,0,0)
        W21ZL=0.
        W21ZN=0.
        W21ZU=0.
        W21ZD=0.
C          ...w1 + f + fbar
      ELSE
        CONST=E**4*(COTW+TANW)**2/96./PI**3/MW2
        UPPER=(MW2**2+MW1**2)/2./MW2
        TMP(1)=MW2
        TMP(2)=MW1
        TMP(3)=AMZ
        TMP(4)=SNW1*SNW2
        TMP(5)=XWINO
        TMP(6)=YWINO
        XINTGL=SSXINT(MW1,SSWWF1,UPPER)
        W21ZL=(APL**2+BTL**2)*CONST*XINTGL
        W21ZN=(APN**2+BTN**2)*CONST*XINTGL
        W21ZU=3*(APU**2+BTU**2)*CONST*XINTGL
        W21ZD=3*(APD**2+BTD**2)*CONST*XINTGL
      END IF
C     do w2 ->w1+q+qbar via sq'
      TMP(1)=MW2
      TMP(3)=MW1
      UPPER=(MW2**2-MW1**2)/2./MW2
      IF (MW2.LT.AMULSS) THEN
        TMP(2)=AMULSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21D=3*(ADWI(2)*ADWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21D=0.
      END IF
      IF (MW2.LT.AMDLSS) THEN
        TMP(2)=AMDLSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21U=3*(AUWI(2)*AUWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21U=0.
      END IF
      IF (MW2.LT.AMCLSS) THEN
        TMP(2)=AMCLSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21S=3*(ADWI(2)*ADWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21S=0.
      END IF
      IF (MW2.LT.AMSLSS) THEN
        TMP(2)=AMSLSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21C=3*(AUWI(2)*AUWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21C=0.
      END IF
      IF (MW2.LT.AMN1SS) THEN
        TMP(2)=AMN1SS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21E=(ALWI(2)*ALWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21E=0.
      END IF
      IF (MW2.LT.AMN2SS) THEN
        TMP(2)=AMN2SS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21M=(ALWI(2)*ALWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21M=0.
      END IF
      IF (MW2.LT.AMN3SS) THEN
        TMP(2)=AMN3SS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21L=(ALWI(2)*ALWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21L=0.
      END IF
      IF (MW2.LT.AMELSS) THEN
        TMP(2)=AMELSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21N1=(ANWI(2)*ANWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21N1=0.
      END IF
      IF (MW2.LT.AMMLSS) THEN
        TMP(2)=AMMLSS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21N2=(ANWI(2)*ANWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21N2=0.
      END IF
C     !!! W2->W1+NU_TAU+NU_TAUBAR NEEDS UPDATING FOR STAU MIXING
      IF (MW2.LT.AML1SS) THEN
        TMP(2)=AML1SS
        XINTGL=SSXINT(0.,SSWZF2,UPPER)
        W21N3=(ANWI(2)*ANWI(1))**2*XINTGL/2./MW2/(2*PI)**5
      ELSE
        W21N3=0.
      END IF
C-----WINO-2 ->WINO-1 +BBBAR NEEDS UPDATE -------------------
C-----WINO-2 ->WINO-1 +TTBAR NEEDS UPDATE -------------------
C-----THESE ALL LACK INTERFERENCE TERMS AS WELL
      W21D=W21D+W21ZD
      W21U=W21U+W21ZU
      W21S=W21S+W21ZD
      W21C=W21C+W21ZU
      W21N1=W21N1+W21ZN
      W21N2=W21N2+W21ZN
      W21N3=W21N3+W21ZN
      W21E=W21E+W21ZL
      W21M=W21M+W21ZL
      W21L=W21L+W21ZL
      IF(W21D.GT.0.) THEN
        CALL SSSAVE(ISW2,W21D,ISW1,IDDN,-IDDN,0,0)
      END IF
      IF(MW2.GT.(MW1+2*AMST+1.)) THEN
        CALL SSSAVE(ISW2,W21S,ISW1,IDST,-IDST,0,0)
      END IF
C      IF (MW2.GT.(MW1+2*AMBT+2.)) THEN
C        CALL SSSAVE(ISW2,W21D,ISW1,IDBT,-IDBT,0,0)
C      END IF
      IF(W21U.GT.0.) THEN
        CALL SSSAVE(ISW2,W21U,ISW1,IDUP,-IDUP,0,0)
      END IF
      IF (MW2.GT.(MW1+2*AMCH+1.)) THEN
        CALL SSSAVE(ISW2,W21C,ISW1,IDCH,-IDCH,0,0)
      ENDIF
C        IF (MW2.GT.(MW1+2*AMTP+2.)) THEN
C          CALL SSSAVE(ISW2,W21U,ISW1,IDTP,-IDTP,0,0)
C        END IF
      IF(W21N1.GT.0.) THEN
        CALL SSSAVE(ISW2,W21N1,ISW1,IDNE,-IDNE,0,0)
      ENDIF
      IF(W21N2.GT.0.) THEN
        CALL SSSAVE(ISW2,W21N2,ISW1,IDNM,-IDNM,0,0)
      ENDIF
      IF(W21N3.GT.0.) THEN
        CALL SSSAVE(ISW2,W21N3,ISW1,IDNT,-IDNT,0,0)
      ENDIF
      IF(W21E.GT.0.) THEN
        CALL SSSAVE(ISW2,W21E,ISW1,IDE,-IDE,0,0)
      ENDIF
      IF(MW2.GT.(MW1+2*AMMU+1.)) THEN
        CALL SSSAVE(ISW2,W21M,ISW1,IDMU,-IDMU,0,0)
      ENDIF
      IF(MW2.GT.(MW1+2*AMTAU+1.)) THEN
        CALL SSSAVE(ISW2,W21L,ISW1,IDTAU,-IDTAU,0,0)
      ENDIF
C
C          w2 --> w1 + higgs
C
C          w2 --> w1 + hl
      IF(MW2.GT.(MW1+AMHL)) THEN
        WID=G*G/16./PI/MW2**3*SQRT(SSXLAM(MW2**2,MW1**2,AMHL**2))*
     $  ((SL*SL+PL*PL)*(MW2*MW2+MW1*MW1-AMHL*AMHL)/2.+
     $  (SL*SL-PL*PL)*MW2*MW1)
        CALL SSSAVE(ISW2,WID,ISW1,ISHL,0,0,0)
      ENDIF
C          w2 --> w1 + hh
      IF(MW2.GT.(MW1+AMHH)) THEN
        WID=G*G/16./PI/MW2**3*SQRT(SSXLAM(MW2**2,MW1**2,AMHH**2))*
     $  ((SH*SH+PH*PH)*(MW2*MW2+MW1*MW1-AMHH*AMHH)/2.+
     $  (SH*SH-PH*PH)*MW2*MW1)
        CALL SSSAVE(ISW2,WID,ISW1,ISHH,0,0,0)
      ENDIF
C          w2 --> w1 + ha
      IF(MW2.GT.(MW1+AMHA)) THEN
        WID=G*G/16./PI/MW2**3*SQRT(SSXLAM(MW2**2,MW1**2,AMHA**2))*
     $  ((SP*SP+PP*PP)*(MW2*MW2+MW1*MW1-AMHA*AMHA)/2.+
     $  (SP*SP-PP*PP)*MW2*MW1)
        CALL SSSAVE(ISW2,WID,ISW1,ISHA,0,0,0)
      END IF
C
C          Normalize wi branching ratios
C
      CALL SSNORM(ISW1)
      CALL SSNORM(ISW2)
C
      RETURN
      END
