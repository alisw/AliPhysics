#include "isajet/pilot.h"
        SUBROUTINE SSZIBF
C-----------------------------------------------------------------------
C       This subroutine calculates the neutralino (zi) subset of
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
      INTEGER ISZ(4),THIZ,THJZ
C
      REAL MWIW,SL,PP,SP,PL,MZIZ,PH,SH,COSB,SINB,COSA,SINA
     $,MZJZ,FAC,COSBE,SINBE
      REAL STHW,CTHW
      REAL EF,A,B,TANB,FB,FT,SR2,G,PI,GP,FL
     $,MW2,SNW1,MW1,YM,BE,SNW2,XM,THX,THY
     $,BTN,APD,APL,APU,BTD,BTL,APN,BTU
      REAL TANW,COTW,XWINO,YWINO,SN,SNIW,SNJZ
      REAL SSXINT,SSXLAM
      REAL WID,WIDEE,E,
     $TERMLL,TERMRR,TERMZZ,TERMLZ,TERMRZ,TERMLH
      REAL FACT,ALIZ1,ALJZ1,ALIZ2,ALJZ2,BEIZ1,BEJZ1,BEIZ2,BEJZ2,
     $SGNIJ,XUPPER,GLLF1,GRRF1,GLRF1,GF1,GLLF2,GRRF2,GLRF2,GF2,
     $GLGL,GRGR,GLGR,GRGL,GF12,GF,TERMHL,TERMHH,TERMHA,
     $TERM1Z,TERM2Z,TERM1L,TERM2L,TERM1H,TERM2H,TERM1A,TERM2A
      REAL COSL,SINL,BPWI(2),BPLWI(2)
      REAL BWI(2),AS,BS,COST,SINT,SNIZ
      REAL GLLF1L,GRRF1R,GLGLL,GRGRR,GLLF2L,GRRF2R,TMZZRL,TMZZLR,
     $TM1ZRL,TM1ZLR,TM2ZRL,TM2ZLR
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      REAL SUALFE,MTAMTA,MTAMB,MTAMZ,AMPL
      REAL FUDGE
      DOUBLE PRECISION SSMQCD
      COMPLEX ZI,ZONE,ZA,ZB,ZPP,ZPM,ZAUIZ,ZBUIZ,ZADIZ,ZBDIZ,
     $ZALIZ,ZBLIZ,Z1(2),Z2(2)
      INTEGER IW,JZ,IZ,ISZIZ,ISWIW
      DATA FUDGE/1.0/
      DATA ZONE/(1.,0.)/,ZI/(0.,1.)/
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
C-----------------------------------------------------------------------
C          Generate Neutralino zi Branching Fractions
C-----------------------------------------------------------------------
      ISZ(1)=ISZ1
      ISZ(2)=ISZ2
      ISZ(3)=ISZ3
      ISZ(4)=ISZ4
C
      DO 200 IZ=2,4
        MZIZ=ABS(AMZISS(IZ))
        SNIZ=SIGN(1.,AMZISS(IZ))
        IF (SNIZ.EQ.1.) THEN
           THIZ=0
        ELSE
           THIZ=1
        END IF
        ISZIZ=ISZ(IZ)
C          z2 --> z1 + photon
C          !!! NEEDS UPDATING
        IF (IZ.EQ.2.AND.(ABS(ZMIXSS(3,1)).LE.1.E-4).AND.
     $  (ABS(ZMIXSS(4,1)).LE.1.E-4)) THEN
          CALL SSSAVE(ISZ2,1.,ISZ1,IDGM,0,0,0)
          GOTO 200
        END IF
C
C          zi --> wi + x partial widths
C
        DO 205 IW=1,2
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
C          zi -> wj + w
C
          IF(MZIZ.GT.MWIW+AMW) THEN
            EF=MZIZ**2+MWIW**2-AMW**2+((MZIZ**2-MWIW**2)**2
     $         -AMW**4)/AMW/AMW
            WID=G*G*SQRT(SSXLAM(MZIZ**2,AMW**2,MWIW**2))/32./PI/
     $          MZIZ**3*(2.*EF*(XIM(IZ)**2+YIM(IZ)**2)-12.*
     $          MZIZ*MWIW*(XIM(IZ)**2-YIM(IZ)**2))
            CALL SSSAVE(ISZIZ,WID,+ISWIW,-IDW,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISWIW,+IDW,0,0,0)
C
C          zi -> wj + f + fbar (w forbidden) ONLY W CONTRIBUTION INCLUDED!
C
          ELSEIF(MZIZ.GT.FUDGE*MWIW) THEN
            TMP(1)=MWIW
            TMP(2)=MZIZ
            IF(IW.EQ.1) THEN
              TMP(3)=XIM(IZ)
              TMP(4)=YIM(IZ)
            ELSE
              TMP(3)=XIP(IZ)
              TMP(4)=YIP(IZ)
            ENDIF
            WID=G**4/96./PI**3/MZIZ
     $      *SSXINT(MWIW,SSZWF1,(MWIW**2+MZIZ**2)/2./MZIZ)
            CALL SSSAVE(ISZIZ,3.*WID,-ISWIW,IDUP,-IDDN,0,0)
            Z1(1)=(-ZI)**THIZ*G*TMP(3)
            Z1(2)=(-ZI)**THIZ*G*TMP(4)
            Z2(1)=G/2./SR2
            Z2(2)=-Z2(1)
            CALL SSME3(1,AMW,Z1,Z2)
            CALL SSSAVE(ISZIZ,WID,-ISWIW,IDNE,-IDE,0,0)
            CALL SSME3(1,AMW,Z1,Z2)
            CALL SSSAVE(ISZIZ,WID,-ISWIW,IDNM,-IDMU,0,0)
            CALL SSME3(1,AMW,Z1,Z2)
            CALL SSSAVE(ISZIZ,3.*WID,ISWIW,IDDN,-IDUP,0,0)
            CALL SSME3(1,AMW,Z1,Z2)
            CALL SSSAVE(ISZIZ,WID,ISWIW,IDE,-IDNE,0,0)
            CALL SSME3(1,AMW,Z1,Z2)
            CALL SSSAVE(ISZIZ,WID,ISWIW,IDMU,-IDNM,0,0)
            CALL SSME3(1,AMW,Z1,Z2)
            IF (MZIZ.GT.(MWIW+AMCH+AMST)) THEN
              CALL SSSAVE(ISZIZ,3.*WID,-ISWIW,IDCH,-IDST,0,0)
              CALL SSME3(1,AMW,Z1,Z2)
              CALL SSSAVE(ISZIZ,WID,-ISWIW,IDNT,-IDTAU,0,0)
              CALL SSME3(1,AMW,Z1,Z2)
              CALL SSSAVE(ISZIZ,3.*WID,ISWIW,IDST,-IDCH,0,0)
              CALL SSME3(1,AMW,Z1,Z2)
              CALL SSSAVE(ISZIZ,WID,ISWIW,IDTAU,-IDNT,0,0)
              CALL SSME3(1,AMW,Z1,Z2)
            ENDIF
          ENDIF
C
C          zi --> wj + hc
C
          IF (MZIZ.GT.(MWIW+AMHC)) THEN
          IF (IW.EQ.1) THEN 
            A=(SNW1*COSBE*V2I(IZ)-SIGN(1.,AMZISS(IZ))
     $      *SINBE*V4I(IZ))/2.
            B=(SNW1*COSBE*V2I(IZ)+SIGN(1.,AMZISS(IZ))
     $      *SINBE*V4I(IZ))/2.
          ELSE
            A=(THY*SNW2*COSBE*V1I(IZ)-SIGN(1.,AMZISS(IZ))
     $      *THX*SINBE*V3I(IZ))/2.
            B=(THY*SNW2*COSBE*V1I(IZ)+SIGN(1.,AMZISS(IZ))
     $      *THX*SINBE*V3I(IZ))/2.
          END IF
            WID=SQRT(MWIW**4+MZIZ**4+AMHC**4-2.*(MWIW*MZIZ)**2
     $      -2*(MWIW*AMHC)**2-2*(MZIZ*AMHC)**2)/8./PI/MZIZ**3
     $      *((A*A+B*B)*(MWIW*MWIW+MZIZ*MZIZ-AMHC*AMHC)/2.
     $      +(A*A-B*B)*MWIW*MZIZ)
            CALL SSSAVE(ISZIZ,WID,+ISWIW,-ISHC,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISWIW,+ISHC,0,0,0)
          ENDIF
205     CONTINUE
C
C          zi --> zj + z
C
C          Note that if m(zi) > m(zj) + m(z), then the z terms are
C          omitted from the zi -> zj + f + fbar calculation, so there
C          is no double counting with zi -> zj + z.
C
        DO 210 JZ=1,IZ-1
          MZJZ=ABS(AMZISS(JZ))
          IF(MZIZ.GT.(MZJZ+AMZ)) THEN
            SN=-1.*SIGN(1.,AMZISS(IZ))*SIGN(1.,AMZISS(JZ))
            WID=WIJ(JZ,IZ)**2/(2.*PI)/(MZIZ)**3*SQRT(MZIZ**4+MZJZ**4
     $      +AMZ**4-2.*(MZIZ*MZJZ)**2-2.*(MZIZ*AMZ)**2-2.*(MZJZ*AMZ)**2)
     $      *((MZIZ**2+MZJZ**2-AMZ**2)/2.+((MZIZ**2-MZJZ**2)**2-AMZ**4)/
     $      2./AMZ**2-3.*SN*MZIZ*MZJZ)
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDZ,0,0,0)
          END IF
210     CONTINUE
C
C          zi --> zj + f + fbar
C
        DO 220 JZ=1,IZ-1
          MZJZ=ABS(AMZISS(JZ))
          SNJZ=SIGN(1.,AMZISS(JZ))
          IF (SNJZ.EQ.1.) THEN
           THJZ=0
          ELSE
           THJZ=1
          END IF
          IF(MZIZ.LT.FUDGE*MZJZ) GO TO 220
          FAC=1./2./MZIZ/(2.*PI)**5*PI**2*MZIZ**2
C          Leptons -- Z decay allowed, so omit Z
            TMP(1)=MZIZ
            TMP(2)=MZJZ
            TMP(4)=-SNIJ(JZ,IZ)
            TMP(6)=0.
C          zi -> zj + u + ubar
          IF (MZIZ.GT.(MZJZ+2*AMUP)) THEN
            IF (MZIZ.LT.AMULSS) THEN
            TMP(3)=AMULSS
            TERMLL=3*2*AUI(JZ)**2*AUI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMURSS) THEN
            TMP(3)=AMURSS
            TERMRR=3*2*BUI(JZ)**2*BUI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=3*(APU**2+BTU**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMULSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMULSS
            TERMLZ=3*8*E*(APU-BTU)*WIJ(JZ,IZ)*AUI(JZ)*AUI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMURSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMURSS
            TERMRZ=-3*8*E*(APU+BTU)*WIJ(JZ,IZ)*BUI(JZ)*BUI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDUP,-IDUP,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APU
            Z2(2)=BTU
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*AUI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMULSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BUI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BUI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMURSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*AUI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMULSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BUI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BUI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMURSS,Z1,Z2)
          END IF
C          zi -> zj + d + dbar
          TMP(4)=-SNIJ(JZ,IZ)
          IF (MZIZ.GT.(MZJZ+2*AMDN)) THEN
            IF (MZIZ.LT.AMDLSS) THEN
            TMP(3)=AMDLSS
            TERMLL=3*2*ADI(JZ)**2*ADI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMDRSS) THEN
            TMP(3)=AMDRSS
            TERMRR=3*2*BDI(JZ)**2*BDI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=3*(APD**2+BTD**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMDLSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMDLSS
            TERMLZ=3*8*E*(APD-BTD)*WIJ(JZ,IZ)*ADI(JZ)*ADI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMDRSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMDRSS
            TERMRZ=-3*8*E*(APD+BTD)*WIJ(JZ,IZ)*BDI(JZ)*BDI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDDN,-IDDN,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APD
            Z2(2)=BTD
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ADI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMDLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BDI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BDI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMDRSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ADI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMDLSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BDI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BDI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMDRSS,Z1,Z2)
          END IF
C          zi -> zj + s + sbar
          TMP(4)=-SNIJ(JZ,IZ)
          IF (MZIZ.GT.(MZJZ+2*AMST)) THEN
            IF (MZIZ.LT.AMSLSS) THEN
            TMP(3)=AMSLSS
            TERMLL=3*2*ADI(JZ)**2*ADI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMSRSS) THEN
            TMP(3)=AMSRSS
            TERMRR=3*2*BDI(JZ)**2*BDI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=3*(APD**2+BTD**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMSLSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMSLSS
            TERMLZ=3*8*E*(APD-BTD)*WIJ(JZ,IZ)*ADI(JZ)*ADI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMSRSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMSRSS
            TERMRZ=-3*8*E*(APD+BTD)*WIJ(JZ,IZ)*BDI(JZ)*BDI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDST,-IDST,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APD
            Z2(2)=BTD
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ADI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMSLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BDI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BDI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMSRSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ADI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ADI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMSLSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BDI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BDI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMSRSS,Z1,Z2)
          END IF
C          zi -> zj + c + cbar
          TMP(4)=-SNIJ(JZ,IZ)
          IF (MZIZ.GT.(MZJZ+2*AMCH)) THEN
            IF (MZIZ.LT.AMCLSS) THEN
            TMP(3)=AMCLSS
            TERMLL=3*2*AUI(JZ)**2*AUI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMCRSS) THEN
            TMP(3)=AMCRSS
            TERMRR=3*2*BUI(JZ)**2*BUI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=3*(APU**2+BTU**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMCLSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMCLSS
            TERMLZ=3*8*E*(APU-BTU)*WIJ(JZ,IZ)*AUI(JZ)*AUI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMCRSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMCRSS
            TERMRZ=-3*8*E*(APU+BTU)*WIJ(JZ,IZ)*BUI(JZ)*BUI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDCH,-IDCH,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APU
            Z2(2)=BTU
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*AUI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMCLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BUI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BUI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMCRSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*AUI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*AUI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMCLSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BUI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BUI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMCRSS,Z1,Z2)
          END IF
C          zi -> zj + b+ bbar ; mixing/yukawa effects now included!
C                              thanks to M. Drees
          FACT=1./2./(2*PI)**5/2./MZIZ
          ALIZ1=ADI(IZ)*COSB-FB*ZMIXSS(2,IZ)*SINB
          ALJZ1=ADI(JZ)*COSB-FB*ZMIXSS(2,JZ)*SINB
          ALIZ2=ADI(IZ)*SINB+FB*ZMIXSS(2,IZ)*COSB
          ALJZ2=ADI(JZ)*SINB+FB*ZMIXSS(2,JZ)*COSB
          BEIZ1=BDI(IZ)*SINB+FB*ZMIXSS(2,IZ)*COSB
          BEJZ1=BDI(JZ)*SINB+FB*ZMIXSS(2,JZ)*COSB
          BEIZ2=-BDI(IZ)*COSB+FB*ZMIXSS(2,IZ)*SINB
          BEJZ2=-BDI(JZ)*COSB+FB*ZMIXSS(2,JZ)*SINB
          SGNIJ=-SNIJ(JZ,IZ)
          XUPPER=(MZIZ**2+AMBT**2-(AMBT+MZJZ)**2)/2./MZIZ
          IF (MZIZ.GT.(MZJZ+2*AMBT)) THEN
            TMP(1)=MZIZ
            TMP(2)=AMBT
            TMP(3)=MZJZ
            TMP(4)=AMB1SS
            TMP(5)=AMB1SS
            TMP(6)=AMBT
            IF (MZIZ.LT.AMB1SS) THEN
              GLLF1=4*ALIZ1**2*((ALJZ1**2+BEJZ1**2)*
     ,             SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*ALJZ1**2*
     ,             SSXINT(AMBT,SSGX2,XUPPER))
              GRRF1=4*BEIZ1**2*((ALJZ1**2+BEJZ1**2)*
     ,             SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*BEJZ1**2*
     ,             SSXINT(AMBT,SSGX2,XUPPER))
              GLRF1=-8*ALIZ1*BEIZ1*ALJZ1*BEJZ1*
     ,               SSXINT(AMBT,SSGX8,XUPPER)
              GF1=GLLF1+GRRF1+GLRF1
            ELSE
              GF1=0.
            END IF
            IF (MZIZ.LT.AMB1SS) THEN
            TMP(4)=AMB1SS
            TMP(5)=AMB2SS
              GLGL=8*ALIZ1*ALIZ2*(ALJZ1*ALJZ2+BEJZ1*BEJZ2)*
     ,             SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*8*ALIZ1*ALIZ2*
     ,             ALJZ1*ALJZ2*SSXINT(AMBT,SSGX2,XUPPER)
              GRGR=8*BEIZ1*BEIZ2*(ALJZ1*ALJZ2+BEJZ1*BEJZ2)*
     ,             SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*8*BEIZ1*BEIZ2*
     ,             BEJZ1*BEJZ2*SSXINT(AMBT,SSGX2,XUPPER)
              GLGR=-8*ALIZ1*BEIZ2*ALJZ2*BEJZ1*
     ,              SSXINT(AMBT,SSGX8,XUPPER)
              GRGL=-8*ALIZ2*BEIZ1*ALJZ1*BEJZ2*
     ,              SSXINT(AMBT,SSGX8,XUPPER)
              GF12=GLGL+GRGR+GLGR+GRGL
            ELSE
              GF12=0.
            END IF
            IF (MZIZ.LT.AMB2SS) THEN
            TMP(4)=AMB2SS
            TMP(5)=AMB2SS
              GLLF2=4*ALIZ2**2*((ALJZ2**2+BEJZ2**2)*
     ,              SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*ALJZ2**2*
     ,              SSXINT(AMBT,SSGX2,XUPPER))
              GRRF2=4*BEIZ2**2*((ALJZ2**2+BEJZ2**2)*
     ,              SSXINT(AMBT,SSGX1,XUPPER)+SGNIJ*BEJZ2**2*
     ,              SSXINT(AMBT,SSGX2,XUPPER))
              GLRF2=-8*ALIZ2*BEIZ2*ALJZ2*BEJZ2*
     ,              SSXINT(AMBT,SSGX8,XUPPER)
              GF2=GLLF2+GRRF2+GLRF2
            ELSE
              GF2=0.
            END IF
            GF=FACT*(GF1+GF2+GF12)
            TMP(2)=MZJZ
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APD**2+BTD**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $       *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2-4*AMBT**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHL)) THEN
            TMP(3)=AMHL
            TMP(5)=AMHL
            TMP(4)=-SNIJ(JZ,IZ)
            TERMHL=G**2/64./PI**3/MZIZ*(MBQ*SINA*(XLIJ(JZ,IZ)+
     $        XLIJ(IZ,JZ))/AMW/COSBE)**2*
     $        SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMBT**2)/2./MZIZ)
            ELSE
            TERMHL=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHH)) THEN
            TMP(3)=AMHH
            TMP(5)=AMHH
            TMP(4)=-SNIJ(JZ,IZ)
            TERMHH=G**2/64./PI**3/MZIZ*(MBQ*COSA*(XHIJ(JZ,IZ)+
     $       XHIJ(IZ,JZ))/AMW/COSBE)**2*
     $       SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMBT**2)/2./MZIZ)
            ELSE
            TERMHH=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHH).AND.MZIZ.LE.(MZJZ+AMHL)) THEN
            TMP(3)=AMHL
            TMP(5)=AMHH
            TMP(4)=-SNIJ(JZ,IZ)
            TERMLH=2*G**2/64./PI**3/MZIZ*(MBQ/AMW/COSBE)**2*
     $            (COSA*(XHIJ(JZ,IZ)+XHIJ(IZ,JZ))*SINA*
     $       (XLIJ(JZ,IZ)+XLIJ(IZ,JZ)))*
     $       SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMBT**2)/2./MZIZ)
            ELSE
            TERMLH=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHA)) THEN
            TMP(3)=AMHA
            TMP(5)=AMHA
            TMP(4)=SNIJ(JZ,IZ)
            TERMHA=G**2*TANB**2/64./PI**3/MZIZ*(MBQ*(XPIJ(JZ,IZ)+
     $            XPIJ(IZ,JZ))/AMW)**2*
     $       SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMBT**2)/2./MZIZ)
            ELSE
            TERMHA=0.
            END IF
            IF (MZIZ.LT.AMB1SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=-SGNIJ
            TMP(3)=AMB1SS
            TERM1Z=8*E*WIJ(JZ,IZ)*(ALJZ1*ALIZ1*(APD-BTD)-
     $       BEJZ1*BEIZ1*(APD+BTD))/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMBT**2,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERM1Z=0.
            END IF
            IF (MZIZ.LT.AMB2SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=-SGNIJ
            TMP(3)=AMB2SS
            TERM2Z=8*E*WIJ(JZ,IZ)*(ALJZ2*ALIZ2*(APD-BTD)-
     $       BEJZ2*BEIZ2*(APD+BTD))/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMBT**2,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERM2Z=0.
            END IF
            IF (MZIZ.LT.AMB1SS.AND.MZIZ.LT.(MZJZ+AMHL)) THEN
            TMP(3)=AMB1SS
            TMP(4)=SGNIJ
            TMP(5)=AMHL
            TERM1L=2*PI*PI*G*MBQ*SINA*SGNIJ*(XLIJ(IZ,JZ)+
     $        XLIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $       /MZIZ/AMW/COSBE*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1L=FACT*TERM1L
            ELSE
            TERM1L=0.
            END IF
            IF (MZIZ.LT.AMB2SS.AND.MZIZ.LT.(MZJZ+AMHL)) THEN
            TMP(3)=AMB2SS
            TMP(4)=SGNIJ
            TMP(5)=AMHL
            TERM2L=2*PI*PI*G*MBQ*SINA*SGNIJ*(XLIJ(IZ,JZ)+
     $        XLIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $       /MZIZ/AMW/COSBE*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2L=FACT*TERM2L
            ELSE
            TERM2L=0.
            END IF
            IF (MZIZ.LT.AMB1SS.AND.MZIZ.LT.(MZJZ+AMHH)) THEN
            TMP(3)=AMB1SS
            TMP(4)=SGNIJ
            TMP(5)=AMHH
            TERM1H=2*PI*PI*G*MBQ*COSA*SGNIJ*(XHIJ(IZ,JZ)+
     $        XHIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $       /MZIZ/AMW/COSBE*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1H=FACT*TERM1H
            ELSE
            TERM1H=0.
            END IF
            IF (MZIZ.LT.AMB2SS.AND.MZIZ.LT.(MZJZ+AMHH)) THEN
            TMP(3)=AMB2SS
            TMP(4)=SGNIJ
            TMP(5)=AMHH
            TERM2H=2*PI*PI*G*MBQ*COSA*SGNIJ*(XHIJ(IZ,JZ)+
     $        XHIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $       /MZIZ/AMW/COSBE*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2H=FACT*TERM2H
            ELSE
            TERM2H=0.
            END IF
            IF (MZIZ.LT.AMB1SS.AND.MZIZ.LT.(MZJZ+AMHA)) THEN
            TMP(3)=AMB1SS
            TMP(4)=-SGNIJ
            TMP(5)=AMHA
            TERM1A=-2*PI*PI*G*MBQ*TANB*SGNIJ*(XPIJ(IZ,JZ)+
     $        XPIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $        /MZIZ/AMW*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1A=FACT*TERM1A
            ELSE
            TERM1A=0.
            END IF
            IF (MZIZ.LT.AMB2SS.AND.MZIZ.LT.(MZJZ+AMHA)) THEN
            TMP(3)=AMB2SS
            TMP(4)=-SGNIJ
            TMP(5)=AMHA
            TERM2A=-2*PI*PI*G*MBQ*TANB*SGNIJ*(XPIJ(IZ,JZ)+
     $        XPIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $        /MZIZ/AMW*SSXINT(4*AMBT**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2A=FACT*TERM2A
            ELSE
            TERM2A=0.
            END IF
            WID=3*(GF+TERMZZ+TERMHL+TERMHH+TERMLH+TERMHA+TERM1Z+
     $        TERM2Z+TERM1L+TERM2L+TERM1H+TERM2H+TERM1A+TERM2A)
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDBT,-IDBT,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APD
            Z2(2)=BTD
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=((ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-(ZI**THIZ*BDI(IZ)-
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
            Z1(2)=((-ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-(ZI**THIZ*BDI(IZ)+
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.
            Z2(1)=CONJG(((ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*COSB-(ZI**THJZ*BDI(JZ)-
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINB)/2.)
            Z2(2)=-CONJG(((-ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*COSB-(ZI**THJZ*BDI(JZ)+
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINB)/2.)
            CALL SSME3(2,AMB1SS,Z1,Z2)
            Z1(1)=((ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+(ZI**THIZ*BDI(IZ)-
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
            Z1(2)=((-ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+(ZI**THIZ*BDI(IZ)+
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.
            Z2(1)=CONJG(((ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*SINB+(ZI**THJZ*BDI(JZ)-
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSB)/2.)
            Z2(2)=-CONJG(((-ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*SINB+(ZI**THJZ*BDI(JZ)+
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSB)/2.)
            CALL SSME3(2,AMB2SS,Z1,Z2)
            Z1(1)=((ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*COSB-(ZI**THJZ*BDI(JZ)-
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINB)/2.
            Z1(2)=((-ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*COSB-(ZI**THJZ*BDI(JZ)+
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINB)/2.
            Z2(1)=-CONJG(((ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-(ZI**THIZ*BDI(IZ)-
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.)
            Z2(2)=CONJG(((-ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB-(ZI**THIZ*BDI(IZ)+
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB)/2.)
            CALL SSME3(3,AMB1SS,Z1,Z2)
            Z1(1)=((ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*SINB+(ZI**THJZ*BDI(JZ)-
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSB)/2.
            Z1(2)=((-ZI**(THJZ)*(-1.)**(THJZ+1)*ADI(JZ)-
     $            FB*ZMIXSS(2,JZ)*ZI**THJZ)*SINB+(ZI**THJZ*BDI(JZ)+
     $            FB*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSB)/2.
            Z2(1)=-CONJG(((ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+(ZI**THIZ*BDI(IZ)-
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.)
            Z2(2)=CONJG(((-ZI**(THIZ)*(-1.)**(THIZ+1)*ADI(IZ)-
     $            FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB+(ZI**THIZ*BDI(IZ)+
     $            FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB)/2.)
            CALL SSME3(3,AMB2SS,Z1,Z2)
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=(-ZI)**(THIZ+THJZ)*(XLIJ(IZ,JZ)+XLIJ(JZ,IZ))
             Z1(2)=0.
            ELSE
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ)*(XLIJ(IZ,JZ)+XLIJ(JZ,IZ))
            END IF
            Z2(1)=-G*MBQ*SINA/2./AMW/COSBE
            Z2(2)=0.
            CALL SSME3(4,AMHL,Z1,Z2)
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=(-ZI)**(THIZ+THJZ)*(XHIJ(IZ,JZ)+XHIJ(JZ,IZ))
             Z1(2)=0.
            ELSE
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ)*(XHIJ(IZ,JZ)+XHIJ(JZ,IZ))
            END IF
            Z2(1)=-G*MBQ*COSA/2./AMW/COSBE
            Z2(2)=0.
            CALL SSME3(4,AMHH,Z1,Z2)
            IF (-SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ+1)*(XPIJ(IZ,JZ)+XPIJ(JZ,IZ))
            ELSE
             Z1(1)=(-ZI)**(THIZ+THJZ+1)*(XPIJ(IZ,JZ)+XPIJ(JZ,IZ))
             Z1(2)=0.
            END IF
            Z2(1)=0.
            Z2(2)=ZI*G*MBQ*TANB/2./AMW
            CALL SSME3(4,AMHA,Z1,Z2)
          END IF
C          zi -> zj + e + ebar
          TMP(2)=MZJZ
          TMP(4)=-SNIJ(JZ,IZ)
          TMP(6)=0.
          IF (MZIZ.GT.(MZJZ+2*AME)) THEN
            IF (MZIZ.LT.AMELSS) THEN
            TMP(3)=AMELSS
            TERMLL=2*ALI(JZ)**2*ALI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMERSS) THEN
            TMP(3)=AMERSS
            TERMRR=2*BLI(JZ)**2*BLI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APL**2+BTL**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMELSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMELSS
            TERMLZ=8*E*(APL-BTL)*WIJ(JZ,IZ)*ALI(JZ)*ALI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMERSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMERSS
            TERMRZ=-8*E*(APL+BTL)*WIJ(JZ,IZ)*BLI(JZ)*BLI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDE,-IDE,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APL
            Z2(2)=BTL
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ALI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMELSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BLI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BLI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMERSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ALI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMELSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BLI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BLI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMERSS,Z1,Z2)
          END IF
C          zi -> zj + mu + mubar
          TMP(4)=-SNIJ(JZ,IZ)
          IF (MZIZ.GT.(MZJZ+2*AMMU)) THEN
            IF (MZIZ.LT.AMMLSS) THEN
            TMP(3)=AMMLSS
            TERMLL=2*ALI(JZ)**2*ALI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.AMMRSS) THEN
            TMP(3)=AMMRSS
            TERMRR=2*BLI(JZ)**2*BLI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMRR=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APL**2+BTL**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMMLSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMMLSS
            TERMLZ=8*E*(APL-BTL)*WIJ(JZ,IZ)*ALI(JZ)*ALI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            IF (MZIZ.LT.AMMRSS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMMRSS
            TERMRZ=-8*E*(APL+BTL)*WIJ(JZ,IZ)*BLI(JZ)*BLI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMRZ=0.
            END IF
            WID=TERMLL+TERMRR+TERMZZ+TERMLZ+TERMRZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDMU,-IDMU,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APL
            Z2(2)=BTL
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ALI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMMLSS,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*BLI(IZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*BLI(JZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(2,AMMRSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ALI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ALI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMMLSS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*BLI(JZ)/2.
            Z1(2)=Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*BLI(IZ)/2.)
            Z2(2)=-Z2(1)
            CALL SSME3(3,AMMRSS,Z1,Z2)
          END IF
C          zi -> zj + tau + taubar.
C          Mixing/yukawa effects now included thanks to M. Drees
          ALIZ1=-ALI(IZ)*COSL-FL*ZMIXSS(2,IZ)*SINL
          ALJZ1=-ALI(JZ)*COSL-FL*ZMIXSS(2,JZ)*SINL
          ALIZ2=-ALI(IZ)*SINL+FL*ZMIXSS(2,IZ)*COSL
          ALJZ2=-ALI(JZ)*SINL+FL*ZMIXSS(2,JZ)*COSL
          BEIZ1=BLI(IZ)*SINL+FL*ZMIXSS(2,IZ)*COSL
          BEJZ1=BLI(JZ)*SINL+FL*ZMIXSS(2,JZ)*COSL
          BEIZ2=-BLI(IZ)*COSL+FL*ZMIXSS(2,IZ)*SINL
          BEJZ2=-BLI(JZ)*COSL+FL*ZMIXSS(2,JZ)*SINL
          SGNIJ=-SNIJ(JZ,IZ)
          XUPPER=(MZIZ**2+AMTAU**2-(AMTAU+MZJZ)**2)/2./MZIZ
C          Polarization for stau_i -> z2ss+tau, z3ss+tau, z4ss+tau.
C          See below for z1ss+tau.
          IF(JZ.EQ.1) THEN
            PTAU1(IZ)=(BEIZ1**2-ALIZ1**2)/(BEIZ1**2+ALIZ1**2)
            PTAU2(IZ)=(BEIZ2**2-ALIZ2**2)/(BEIZ2**2+ALIZ2**2)
          ENDIF
          IF (MZIZ.GT.(MZJZ+2*AMTAU)) THEN
            TMP(1)=MZIZ
            TMP(2)=AMTAU
            TMP(3)=MZJZ
            TMP(4)=AML1SS
            TMP(5)=AML1SS
            TMP(6)=AMTAU
            IF (MZIZ.LT.AML1SS) THEN
              GLLF1=4*ALIZ1**2*((ALJZ1**2+BEJZ1**2)*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*ALJZ1**2*
     ,             SSXINT(AMTAU,SSGX2,XUPPER))
              GLLF1L=FACT*4*ALIZ1**2*(ALJZ1**2*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*ALJZ1**2*
     ,             SSXINT(AMTAU,SSGX2,XUPPER))
              GRRF1=4*BEIZ1**2*((ALJZ1**2+BEJZ1**2)*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*BEJZ1**2*
     ,             SSXINT(AMTAU,SSGX2,XUPPER))
              GRRF1R=FACT*4*BEIZ1**2*(BEJZ1**2*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*BEJZ1**2*
     ,             SSXINT(AMTAU,SSGX2,XUPPER))
              GLRF1=-8*ALIZ1*BEIZ1*ALJZ1*BEJZ1*
     ,               SSXINT(AMTAU,SSGX8,XUPPER)
              GF1=GLLF1+GRRF1+GLRF1
            ELSE
              GLLF1L=0.
              GRRF1R=0.
              GF1=0.
            END IF
            IF (MZIZ.LT.AML1SS) THEN
            TMP(4)=AML1SS
            TMP(5)=AML2SS
              GLGL=8*ALIZ1*ALIZ2*(ALJZ1*ALJZ2+BEJZ1*BEJZ2)*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*8*ALIZ1*ALIZ2*
     ,             ALJZ1*ALJZ2*SSXINT(AMTAU,SSGX2,XUPPER)
              GLGLL=FACT*(8*ALIZ1*ALIZ2*ALJZ1*ALJZ2*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*8*ALIZ1*ALIZ2*
     ,             ALJZ1*ALJZ2*SSXINT(AMTAU,SSGX2,XUPPER))
              GRGR=8*BEIZ1*BEIZ2*(ALJZ1*ALJZ2+BEJZ1*BEJZ2)*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*8*BEIZ1*BEIZ2*
     ,             BEJZ1*BEJZ2*SSXINT(AMTAU,SSGX2,XUPPER)
              GRGRR=FACT*(8*BEIZ1*BEIZ2*BEJZ1*BEJZ2*
     ,             SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*8*BEIZ1*BEIZ2*
     ,             BEJZ1*BEJZ2*SSXINT(AMTAU,SSGX2,XUPPER))
              GLGR=-8*ALIZ1*BEIZ2*ALJZ2*BEJZ1*
     ,              SSXINT(AMTAU,SSGX8,XUPPER)
              GRGL=-8*ALIZ2*BEIZ1*ALJZ1*BEJZ2*
     ,              SSXINT(AMTAU,SSGX8,XUPPER)
              GF12=GLGL+GRGR+GLGR+GRGL
            ELSE
              GLGLL=0.
              GRGRR=0.
              GF12=0.
            END IF
            IF (MZIZ.LT.AML2SS) THEN
            TMP(4)=AML2SS
            TMP(5)=AML2SS
              GLLF2=4*ALIZ2**2*((ALJZ2**2+BEJZ2**2)*
     ,              SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*ALJZ2**2*
     ,              SSXINT(AMTAU,SSGX2,XUPPER))
              GLLF2L=FACT*4*ALIZ2**2*(ALJZ2**2*
     ,              SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*ALJZ2**2*
     ,              SSXINT(AMTAU,SSGX2,XUPPER))
              GRRF2=4*BEIZ2**2*((ALJZ2**2+BEJZ2**2)*
     ,              SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*BEJZ2**2*
     ,              SSXINT(AMTAU,SSGX2,XUPPER))
              GRRF2R=FACT*4*BEIZ2**2*(BEJZ2**2*
     ,              SSXINT(AMTAU,SSGX1,XUPPER)+SGNIJ*BEJZ2**2*
     ,              SSXINT(AMTAU,SSGX2,XUPPER))
              GLRF2=-8*ALIZ2*BEIZ2*ALJZ2*BEJZ2*
     ,              SSXINT(AMTAU,SSGX8,XUPPER)
              GF2=GLLF2+GRRF2+GLRF2
            ELSE
              GLLF2L=0.
              GRRF2R=0.
              GF2=0.
            END IF
            GF=FACT*(GF1+GF2+GF12)
            TMP(2)=MZJZ
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APL**2+BTL**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $      *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2-4*AMTAU**2)/2./MZIZ)
            TMZZRL=TERMZZ*(APL+BTL)**2/2./(APL**2+BTL**2)
            TMZZLR=TERMZZ*(APL-BTL)**2/2./(APL**2+BTL**2)
            ELSE
            TERMZZ=0.
            TMZZRL=0.
            TMZZLR=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHL)) THEN
            TMP(3)=AMHL
            TMP(5)=AMHL
            TMP(4)=-SNIJ(JZ,IZ)
            TERMHL=G**2/64./PI**3/MZIZ*(MTAMZ*SINA*(XLIJ(JZ,IZ)+
     $         XLIJ(IZ,JZ))/AMW/COSBE)**2*
     $      SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMTAU**2)/2./MZIZ)
            ELSE
            TERMHL=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHH)) THEN
            TMP(3)=AMHH
            TMP(5)=AMHH
            TMP(4)=-SNIJ(JZ,IZ)
            TERMHH=G**2/64./PI**3/MZIZ*(MTAMZ*COSA*(XHIJ(JZ,IZ)+
     $            XHIJ(IZ,JZ))/AMW/COSBE)**2*
     $      SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMTAU**2)/2./MZIZ)
            ELSE
            TERMHH=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHH).AND.MZIZ.LE.(MZJZ+AMHL)) THEN
            TMP(3)=AMHL
            TMP(5)=AMHH
            TMP(4)=-SNIJ(JZ,IZ)
            TERMLH=2*G**2/64./PI**3/MZIZ*(MTAMZ/AMW/COSBE)**2*
     $            (COSA*(XHIJ(JZ,IZ)+XHIJ(IZ,JZ))*SINA*
     $             (XLIJ(JZ,IZ)+XLIJ(IZ,JZ)))*
     $      SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMTAU**2)/2./MZIZ)
            ELSE
            TERMLH=0.
            END IF
            IF (MZIZ.LE.(MZJZ+AMHA)) THEN
            TMP(3)=AMHA
            TMP(5)=AMHA
            TMP(4)=SNIJ(JZ,IZ)
            TERMHA=G**2*TANB**2/64./PI**3/MZIZ*(MTAMZ*(XPIJ(JZ,IZ)+
     $            XPIJ(IZ,JZ))/AMW)**2*
     $      SSXINT(MZJZ,SSZZF4,(MZIZ**2+MZJZ**2-4*AMTAU**2)/2./MZIZ)
            ELSE
            TERMHA=0.
            END IF
            IF (MZIZ.LT.AML1SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=-SGNIJ
            TMP(3)=AML1SS
            TERM1Z=8*E*WIJ(JZ,IZ)*(ALJZ1*ALIZ1*(APL-BTL)-
     $       BEJZ1*BEIZ1*(APL+BTL))/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            TM1ZRL=-8*E*WIJ(JZ,IZ)*BEJZ1*BEIZ1*(APL+BTL)/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            TM1ZLR=8*E*WIJ(JZ,IZ)*ALJZ1*ALIZ1*(APL-BTL)/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERM1Z=0.
            TM1ZRL=0.
            TM1ZLR=0.
            END IF
            IF (MZIZ.LT.AML2SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=-SGNIJ
            TMP(3)=AML2SS
            TERM2Z=8*E*WIJ(JZ,IZ)*(ALJZ2*ALIZ2*(APL-BTL)-
     $       BEJZ2*BEIZ2*(APL+BTL))/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            TM2ZRL=-8*E*WIJ(JZ,IZ)*BEJZ2*BEIZ2*(APL+BTL)/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            TM2ZLR=8*E*WIJ(JZ,IZ)*ALJZ2*ALIZ2*(APL-BTL)/MZIZ
     $      /(2*PI)**5*SSXINT(4*AMTAU**2,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERM2Z=0.
            TM2ZRL=0.
            TM2ZLR=0.
            END IF
            IF (MZIZ.LT.AML1SS.AND.MZIZ.LT.(MZJZ+AMHL)) THEN
            TMP(3)=AML1SS
            TMP(4)=SGNIJ
            TMP(5)=AMHL
            TERM1L=2*PI*PI*G*MTAMZ*SINA*SGNIJ*(XLIJ(IZ,JZ)+
     $        XLIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $      /MZIZ/AMW/COSBE*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1L=FACT*TERM1L
            ELSE
            TERM1L=0.
            END IF
            IF (MZIZ.LT.AML2SS.AND.MZIZ.LT.(MZJZ+AMHL)) THEN
            TMP(3)=AML2SS
            TMP(4)=SGNIJ
            TMP(5)=AMHL
            TERM2L=2*PI*PI*G*MTAMZ*SINA*SGNIJ*(XLIJ(IZ,JZ)+
     $        XLIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $      /MZIZ/AMW/COSBE*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2L=FACT*TERM2L
            ELSE
            TERM2L=0.
            END IF
            IF (MZIZ.LT.AML1SS.AND.MZIZ.LT.(MZJZ+AMHH)) THEN
            TMP(3)=AML1SS
            TMP(4)=SGNIJ
            TMP(5)=AMHH
            TERM1H=2*PI*PI*G*MTAMZ*COSA*SGNIJ*(XHIJ(IZ,JZ)+
     $        XHIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $      /MZIZ/AMW/COSBE*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1H=FACT*TERM1H
            ELSE
            TERM1H=0.
            END IF
            IF (MZIZ.LT.AML2SS.AND.MZIZ.LT.(MZJZ+AMHH)) THEN
            TMP(3)=AML2SS
            TMP(4)=SGNIJ
            TMP(5)=AMHH
            TERM2H=2*PI*PI*G*MTAMZ*COSA*SGNIJ*(XHIJ(IZ,JZ)+
     $        XHIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $    /MZIZ/AMW/COSBE*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2H=FACT*TERM2H
            ELSE
            TERM2H=0.
            END IF
            IF (MZIZ.LT.AML1SS.AND.MZIZ.LT.(MZJZ+AMHA)) THEN
            TMP(3)=AML1SS
            TMP(4)=-SGNIJ
            TMP(5)=AMHA
            TERM1A=-2*PI*PI*G*MTAMZ*TANB*SGNIJ*(XPIJ(IZ,JZ)+
     $        XPIJ(JZ,IZ))*(ALIZ1*BEJZ1+ALJZ1*BEIZ1)
     $       /MZIZ/AMW*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM1A=FACT*TERM1A
            ELSE
            TERM1A=0.
            END IF
            IF (MZIZ.LT.AML2SS.AND.MZIZ.LT.(MZJZ+AMHA)) THEN
            TMP(3)=AML2SS
            TMP(4)=-SGNIJ
            TMP(5)=AMHA
            TERM2A=-2*PI*PI*G*MTAMZ*TANB*SGNIJ*(XPIJ(IZ,JZ)+
     $        XPIJ(JZ,IZ))*(ALIZ2*BEJZ2+ALJZ2*BEIZ2)
     $       /MZIZ/AMW*SSXINT(4*AMTAU**2,SSZZF5,(MZIZ-MZJZ)**2)
            TERM2A=FACT*TERM2A
            ELSE
            TERM2A=0.
            END IF
            WID=GF+TERMZZ+TERMHL+TERMHH+TERMLH+TERMHA+TERM1Z+
     $        TERM2Z+TERM1L+TERM2L+TERM1H+TERM2H+TERM1A+TERM2A
C              tau polarization for 3-body z2 -> z1 tau tau
            IF (IZ.EQ.2.AND.JZ.EQ.1.AND.WID.GT.0.) THEN
              PTAUZZ=(GRRF1R+GRGRR+GRRF2R+TMZZRL+TM1ZRL+TM2ZRL-
     $               (GLLF1L+GLGLL+GLLF2L+TMZZLR+TM1ZLR+TM2ZLR))
     $               /WID
            END IF
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDTAU,-IDTAU,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APL
            Z2(2)=BTL
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=((ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*COSL-(ZI**THIZ*BLI(IZ)-
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINL)/2.
            Z1(2)=((-ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*COSL-(ZI**THIZ*BLI(IZ)+
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINL)/2.
            Z2(1)=CONJG(((ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*COSL-(ZI**THJZ*BLI(JZ)-
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINL)/2.)
            Z2(2)=-CONJG(((-ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*COSL-(ZI**THJZ*BLI(JZ)+
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINL)/2.)
            CALL SSME3(2,AML1SS,Z1,Z2)
            Z1(1)=((ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*SINL+(ZI**THIZ*BLI(IZ)-
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSL)/2.
            Z1(2)=((-ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*SINL+(ZI**THIZ*BLI(IZ)+
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSL)/2.
            Z2(1)=CONJG(((ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*SINL+(ZI**THJZ*BLI(JZ)-
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSL)/2.)
            Z2(2)=-CONJG(((-ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*SINL+(ZI**THJZ*BLI(JZ)+
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSL)/2.)
            CALL SSME3(2,AML2SS,Z1,Z2)
            Z1(1)=((ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*COSL-(ZI**THJZ*BLI(JZ)-
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINL)/2.
            Z1(2)=((-ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*COSL-(ZI**THJZ*BLI(JZ)+
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*SINL)/2.
            Z2(1)=-CONJG(((ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*COSL-(ZI**THIZ*BLI(IZ)-
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINL)/2.)
            Z2(2)=CONJG(((-ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*COSL-(ZI**THIZ*BLI(IZ)+
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINL)/2.)
            CALL SSME3(3,AML1SS,Z1,Z2)
            Z1(1)=((ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*SINL+(ZI**THJZ*BLI(JZ)-
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSL)/2.
            Z1(2)=((-ZI**(THJZ)*(-1.)**(THJZ+1)*ALI(JZ)-
     $            FL*ZMIXSS(2,JZ)*ZI**THJZ)*SINL+(ZI**THJZ*BLI(JZ)+
     $            FL*ZMIXSS(2,JZ)*(-ZI)**THJZ)*COSL)/2.
            Z2(1)=-CONJG(((ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*SINL+(ZI**THIZ*BLI(IZ)-
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSL)/2.)
            Z2(2)=CONJG(((-ZI**(THIZ)*(-1.)**(THIZ+1)*ALI(IZ)-
     $            FL*ZMIXSS(2,IZ)*ZI**THIZ)*SINL+(ZI**THIZ*BLI(IZ)+
     $            FL*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSL)/2.)
            CALL SSME3(3,AML2SS,Z1,Z2)
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=(-ZI)**(THIZ+THJZ)*(XLIJ(IZ,JZ)+XLIJ(JZ,IZ))
             Z1(2)=0.
            ELSE
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ)*(XLIJ(IZ,JZ)+XLIJ(JZ,IZ))
            END IF
            Z2(1)=-G*MTAMZ*SINA/2./AMW/COSBE
            Z2(2)=0.
            CALL SSME3(4,AMHL,Z1,Z2)
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=(-ZI)**(THIZ+THJZ)*(XHIJ(IZ,JZ)+XHIJ(JZ,IZ))
             Z1(2)=0.
            ELSE
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ)*(XHIJ(IZ,JZ)+XHIJ(JZ,IZ))
            END IF
            Z2(1)=-G*MTAMZ*COSA/2./AMW/COSBE
            Z2(2)=0.
            CALL SSME3(4,AMHH,Z1,Z2)
            IF (-SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=(-ZI)**(THIZ+THJZ+1)*(XPIJ(IZ,JZ)+XPIJ(JZ,IZ))
            ELSE
             Z1(1)=(-ZI)**(THIZ+THJZ+1)*(XPIJ(IZ,JZ)+XPIJ(JZ,IZ))
             Z1(2)=0.
            END IF
            Z2(1)=0.
            Z2(2)=ZI*G*MTAMZ*TANB/2./AMW
            CALL SSME3(4,AMHA,Z1,Z2)
          END IF
C          zi -> zj + nu_e + nu_e bar
          TMP(6)=0.
          IF (MZIZ.GT.MZJZ) THEN
            IF (MZIZ.LT.AMN1SS) THEN
            TMP(4)=-SNIJ(JZ,IZ)
            TMP(3)=AMN1SS
            TERMLL=2*ANI(JZ)**2*ANI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APN**2+BTN**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMN1SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMN1SS
            TERMLZ=8*E*(APN-BTN)*WIJ(JZ,IZ)*ANI(JZ)*ANI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            WID=TERMLL+TERMZZ+TERMLZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDNE,-IDNE,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APN
            Z2(2)=BTN
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMN1SS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMN1SS,Z1,Z2)
          END IF
C          zi -> zj + nu_mu + nu_mu bar
          IF (MZIZ.GT.MZJZ) THEN
            IF (MZIZ.LT.AMN2SS) THEN
            TMP(4)=-SNIJ(JZ,IZ)
            TMP(3)=AMN2SS
            TERMLL=2*ANI(JZ)**2*ANI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APN**2+BTN**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMN2SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMN2SS
            TERMLZ=8*E*(APN-BTN)*WIJ(JZ,IZ)*ANI(JZ)*ANI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            WID=TERMLL+TERMZZ+TERMLZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDNM,-IDNM,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APN
            Z2(2)=BTN
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMN2SS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMN2SS,Z1,Z2)
          END IF
C          zi -> zj + nu_tau + nu_tau bar
          IF (MZIZ.GT.MZJZ) THEN
            IF (MZIZ.LT.AMN3SS) THEN
            TMP(4)=-SNIJ(JZ,IZ)
            TMP(3)=AMN3SS
            TERMLL=2*ANI(JZ)**2*ANI(IZ)**2*FAC*SSXINT(0.,SSZZF1,1.)
            ELSE
            TERMLL=0.
            END IF
            IF (MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TERMZZ=(APN**2+BTN**2)*E**2*WIJ(JZ,IZ)**2/MZIZ
     $             *SSXINT(MZJZ,SSZZF2,(MZIZ**2+MZJZ**2)/2./MZIZ)
            ELSE
            TERMZZ=0.
            END IF
            IF (MZIZ.LT.AMN3SS.AND.MZIZ.LT.(MZJZ+AMZ)) THEN
            TMP(4)=+SNIJ(JZ,IZ)
            TMP(3)=AMN3SS
            TERMLZ=8*E*(APN-BTN)*WIJ(JZ,IZ)*ANI(JZ)*ANI(IZ)/MZIZ
     $      /(2*PI)**5*SSXINT(0.,SSZZF3,(MZIZ-MZJZ)**2)
            ELSE
            TERMLZ=0.
            END IF
            WID=TERMLL+TERMZZ+TERMLZ
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),IDNT,-IDNT,0,0)
C     Enter information for decay matrix element
            IF (-1.*SNIJ(JZ,IZ).GT.0.) THEN
             Z1(1)=0.
             Z1(2)=2*E*WIJ(JZ,IZ)
            ELSE
             Z1(1)=2*E*WIJ(JZ,IZ)
             Z1(2)=0.
            END IF
            Z2(1)=APN
            Z2(2)=BTN
            CALL SSME3(1,AMZ,Z1,Z2)
            Z1(1)=ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=CONJG(ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(2,AMN3SS,Z1,Z2)
            Z1(1)=ZI**(THJZ-1)*(-1.)**(THJZ+1)*ANI(JZ)/2.
            Z1(2)=-Z1(1)
            Z2(1)=-CONJG(ZI**(THIZ-1)*(-1.)**(THIZ+1)*ANI(IZ)/2.)
            Z2(2)=Z2(1)
            CALL SSME3(3,AMN3SS,Z1,Z2)
          END IF
220     CONTINUE
C
C          zi --> zj + higgs
C
        DO 230 JZ=1,IZ-1
C          zi --> zj + hl
          MZJZ=ABS(AMZISS(JZ))
          IF (MZIZ.GT.(MZJZ+AMHL)) THEN
            SN=SIGN(1.,AMZISS(JZ))*SIGN(1.,AMZISS(IZ))
            WID=(XLIJ(JZ,IZ)+XLIJ(IZ,JZ))**2/8./PI/(MZIZ)**3
     $      *SQRT(MZIZ**4+MZJZ**4+AMHL**4-2.*(MZIZ*MZJZ)**2
     $      -2.*(MZIZ*AMHL)**2-2.*(MZJZ*AMHL)**2)*((MZIZ**2+MZJZ**2
     $      -AMHL**2)/2.+SN*MZIZ*MZJZ)
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),ISHL,0,0,0)
          END IF
C          zi --> zj + hh
          IF (MZIZ.GT.(MZJZ+AMHH)) THEN
            SN=SIGN(1.,AMZISS(JZ))*SIGN(1.,AMZISS(IZ))
            WID=(HIJ(JZ,IZ)+HIJ(IZ,JZ))**2/8./PI/(MZIZ)**3
     $      *SQRT(MZIZ**4+MZJZ**4+AMHH**4-2.*(MZIZ*MZJZ)**2
     $      -2.*(MZIZ*AMHH)**2-2.*(MZJZ*AMHH)**2)
     $      *((MZIZ**2+MZJZ**2-AMHH**2)/2.+SN*MZIZ*MZJZ)
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),ISHH,0,0,0)
          END IF
C          zi --> zj + ha
          IF (MZIZ.GT.(MZJZ+AMHA)) THEN
            SN=-SIGN(1.,AMZISS(JZ))*SIGN(1.,AMZISS(IZ))
            WID=(XPIJ(IZ,JZ)+XPIJ(JZ,IZ))**2/8./PI/(MZIZ)**3
     $      *SQRT(MZIZ**4+MZJZ**4+AMHA**4-2.*(MZIZ*MZJZ)**2
     $      -2.*(MZIZ*AMHA)**2-2.*(MZJZ*AMHA)**2)*((MZIZ**2+MZJZ**2
     $      -AMHA**2)/2.+SN*MZIZ*MZJZ)
            CALL SSSAVE(ISZIZ,WID,ISZ(JZ),ISHA,0,0,0)
          END IF
230     CONTINUE
200   CONTINUE
C
C          zi --> squark + qbar; enlarge to include Z1 decays
C                 in case of models with light gravitino
C
      DO 245 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNIZ=SIGN(1.,AMZISS(IZ))
        IF (SNIZ.EQ.1.) THEN
           THIZ=0
        ELSE
           THIZ=1
        END IF
        ISZIZ=ISZ(IZ)
C          left squarks
        IF (MZIZ.GT.(AMULSS+AMUP)) THEN
          WID=3*AUI(IZ)**2*(MZIZ**2+AMUP**2-AMULSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMUP**2,AMULSS**2))
          CALL SSSAVE(ISZIZ,WID,ISUPL,-IDUP,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISUPL,IDUP,0,0,0)
        ENDIF
        IF (MZIZ.GT.(AMDLSS+AMDN)) THEN
          WID=3*ADI(IZ)**2*(MZIZ**2+AMDN**2-AMDLSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMDN**2,AMDLSS**2))
          CALL SSSAVE(ISZIZ,WID,ISDNL,-IDDN,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISDNL,IDDN,0,0,0)
        END IF
        IF (MZIZ.GT.(AMSLSS+AMST)) THEN
          WID=3*ADI(IZ)**2*(MZIZ**2+AMST**2-AMSLSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMST**2,AMSLSS**2))
          CALL SSSAVE(ISZIZ,WID,ISSTL,-IDST,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISSTL,IDST,0,0,0)
        END IF
        IF (MZIZ.GT.(AMCLSS+AMCH)) THEN
          WID=3*AUI(IZ)**2*(MZIZ**2+AMCH**2-AMCLSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMCH**2,AMCLSS**2))
          CALL SSSAVE(ISZIZ,WID,ISCHL,-IDCH,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISCHL,IDCH,0,0,0)
        ENDIF
C          right squarks
        IF (MZIZ.GT.(AMURSS+AMUP)) THEN
          WID=3*BUI(IZ)**2*(MZIZ**2+AMUP**2-AMURSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMUP**2,AMURSS**2))
          CALL SSSAVE(ISZIZ,WID,ISUPR,-IDUP,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISUPR,IDUP,0,0,0)
        END IF
        IF (MZIZ.GT.(AMDRSS+AMDN)) THEN
          WID=3*BDI(IZ)**2*(MZIZ**2+AMDN**2-AMDRSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMDN**2,AMDRSS**2))
          CALL SSSAVE(ISZIZ,WID,ISDNR,-IDDN,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISDNR,IDDN,0,0,0)
        END IF
        IF (MZIZ.GT.(AMSRSS+AMST)) THEN
          WID=3*BDI(IZ)**2*(MZIZ**2+AMST**2-AMSRSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMST**2,AMSRSS**2))
          CALL SSSAVE(ISZIZ,WID,ISSTR,-IDST,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISSTR,IDST,0,0,0)
        END IF
        IF(MZIZ.GT.(AMCRSS+AMCH)) THEN
          WID=3*BUI(IZ)**2*(MZIZ**2+AMCH**2-AMCRSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMCH**2,AMCRSS**2))
            CALL SSSAVE(ISZIZ,WID,ISCHR,-IDCH,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISCHR,IDCH,0,0,0)
        ENDIF
C          z_i --> sbottom + bottom
C
        ZAUIZ=ZI**(THIZ-1)*SNIZ
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZADIZ=ZI**(THIZ-1)*SNIZ
     $  *(G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBDIZ=-2*ZI**(THIZ-1)*GP*ZMIXSS(4,IZ)/3./SR2
        ZALIZ=ZI**(THIZ-1)*SNIZ
     $  *(G/SR2*ZMIXSS(3,IZ)+GP/SR2*ZMIXSS(4,IZ))
        ZBLIZ=-1*ZI**(THIZ-1)*SR2*GP*ZMIXSS(4,IZ)
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        IF(MZIZ.GT.(AMB1SS+AMBT)) THEN
          ZA=((ZI*ZADIZ-ZPP*FB*ZMIXSS(2,IZ))*COSB
     $     -(ZI*ZBDIZ-ZPM*FB*ZMIXSS(2,IZ))*SINB)/2.
          ZB=((-ZI*ZADIZ-ZPP*FB*ZMIXSS(2,IZ))*COSB
     $     -(ZI*ZBDIZ+ZPM*FB*ZMIXSS(2,IZ))*SINB)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=3*(AS*((AMBT+MZIZ)**2-AMB1SS**2)+BS*((MZIZ-AMBT)**2-
     $     AMB1SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AMB1SS**2,AMBT**2))
            CALL SSSAVE(ISZIZ,WID,ISBT1,-IDBT,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISBT1,IDBT,0,0,0)
        ENDIF
        IF(MZIZ.GT.(AMB2SS+AMBT)) THEN
          ZA=((ZI*ZADIZ-ZPP*FB*ZMIXSS(2,IZ))*SINB
     $     +(ZI*ZBDIZ-ZPM*FB*ZMIXSS(2,IZ))*COSB)/2.
          ZB=((-ZI*ZADIZ-ZPP*FB*ZMIXSS(2,IZ))*SINB
     $     +(ZI*ZBDIZ+ZPM*FB*ZMIXSS(2,IZ))*COSB)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=3*(AS*((AMBT+MZIZ)**2-AMB2SS**2)+BS*((MZIZ-AMBT)**2-
     $     AMB2SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AMB2SS**2,AMBT**2))
            CALL SSSAVE(ISZIZ,WID,ISBT2,-IDBT,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISBT2,IDBT,0,0,0)
        ENDIF
C          z_i --> stop + top
C
        IF(MZIZ.GT.AMT1SS+AMTP) THEN
          ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $     -(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
          ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $     -(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=3*(AS*((AMTP+MZIZ)**2-AMT1SS**2)+BS*((MZIZ-AMTP)**2-
     $     AMT1SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AMT1SS**2,AMTP**2))
          CALL SSSAVE(ISZIZ,WID,ISTP1,-IDTP,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISTP1,IDTP,0,0,0)
        ENDIF
        IF(MZIZ.GT.AMT2SS+AMTP) THEN
          ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $     +(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
          ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $     +(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=3*(AS*((AMTP+MZIZ)**2-AMT2SS**2)+BS*((MZIZ-AMTP)**2-
     $     AMT2SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AMT2SS**2,AMTP**2))
          CALL SSSAVE(ISZIZ,WID,ISTP2,-IDTP,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISTP2,IDTP,0,0,0)
        ENDIF
C
C          zi --> slepton + lepton
C
        IF(MZIZ.GT.(AMELSS+AME)) THEN
          WID=ALI(IZ)**2*(MZIZ**2+AME**2-AMELSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AME**2,AMELSS**2))
          CALL SSSAVE(ISZIZ,WID,ISEL,-IDE,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISEL,IDE,0,0,0)
        END IF
        IF(MZIZ.GT.(AMMLSS+AMMU)) THEN
          WID=ALI(IZ)**2*(MZIZ**2+AMMU**2-AMMLSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMMU**2,AMMLSS**2))
          CALL SSSAVE(ISZIZ,WID,ISMUL,-IDMU,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISMUL,IDMU,0,0,0)
        END IF
        IF(MZIZ.GT.(AMERSS+AME)) THEN
          WID=BLI(IZ)**2*(MZIZ**2+AME**2-AMERSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AME**2,AMERSS**2))
          CALL SSSAVE(ISZIZ,WID,ISER,-IDE,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISER,IDE,0,0,0)
        END IF
        IF(MZIZ.GT.(AMMRSS+AMMU)) THEN
          WID=BLI(IZ)**2*(MZIZ**2+AMMU**2-AMMRSS**2)/MZIZ**3/
     $         32./PI*SQRT(SSXLAM(MZIZ**2,AMMU**2,AMMRSS**2))
          CALL SSSAVE(ISZIZ,WID,ISMUR,-IDMU,0,0,0)
          CALL SSSAVE(ISZIZ,WID,-ISMUR,IDMU,0,0,0)
        END IF
        IF(MZIZ.GT.(AML1SS+AMTAU)) THEN
          ZA=((ZI*ZALIZ-ZPP*FL*ZMIXSS(2,IZ))*COSL
     $     -(ZI*ZBLIZ-ZPM*FL*ZMIXSS(2,IZ))*SINL)/2.
          ZB=((-ZI*ZALIZ-ZPP*FL*ZMIXSS(2,IZ))*COSL
     $     -(ZI*ZBLIZ+ZPM*FL*ZMIXSS(2,IZ))*SINL)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=(AS*((AMTAU+MZIZ)**2-AML1SS**2)+BS*((MZIZ-AMTAU)**2-
     $     AML1SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AML1SS**2,AMTAU**2))
            CALL SSSAVE(ISZIZ,WID,ISTAU1,-IDTAU,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISTAU1,IDTAU,0,0,0)
        END IF
        IF(MZIZ.GT.(AML2SS+AMTAU)) THEN
          ZA=((ZI*ZALIZ-ZPP*FL*ZMIXSS(2,IZ))*SINL
     $     +(ZI*ZBLIZ-ZPM*FL*ZMIXSS(2,IZ))*COSL)/2.
          ZB=((-ZI*ZALIZ-ZPP*FL*ZMIXSS(2,IZ))*SINL
     $     +(ZI*ZBLIZ+ZPM*FL*ZMIXSS(2,IZ))*COSL)/2.
          AS=ZA*CONJG(ZA)
          BS=ZB*CONJG(ZB)
          WID=(AS*((AMTAU+MZIZ)**2-AML2SS**2)+BS*((MZIZ-AMTAU)**2-
     $     AML2SS**2))/16./PI/MZIZ**3*
     $     SQRT(SSXLAM(MZIZ**2,AML2SS**2,AMTAU**2))
            CALL SSSAVE(ISZIZ,WID,ISTAU2,-IDTAU,0,0,0)
            CALL SSSAVE(ISZIZ,WID,-ISTAU2,IDTAU,0,0,0)
        ENDIF
        IF(MZIZ.GT.AMN1SS) THEN
          WID=(MZIZ**2-AMN1SS**2)**2/32./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,ISNEL,-IDNE,0,0,0)
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,-ISNEL,IDNE,0,0,0)
        END IF
        IF(MZIZ.GT.AMN2SS) THEN
          WID=(MZIZ**2-AMN2SS**2)**2/32./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,ISNML,-IDNM,0,0,0)
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,-ISNML,IDNM,0,0,0)
        END IF
        IF(MZIZ.GT.AMN3SS) THEN
          WID=(MZIZ**2-AMN3SS**2)**2/32./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,ISNTL,-IDNT,0,0,0)
          CALL SSSAVE(ISZIZ,ANI(IZ)**2*WID,-ISNTL,IDNT,0,0,0)
        END IF
245   CONTINUE
C
C          Compute decays to gravitino for GMSB models
C
      DO 250 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        ISZIZ=ISZ(IZ)
        IF (MZIZ.GT.AMGVSS) THEN
        WID=(ZMIXSS(4,IZ)*CTHW+ZMIXSS(3,IZ)*STHW)**2*MZIZ**5/
     ,       48./PI/(AMGVSS*AMPL)**2
        CALL SSSAVE(ISZIZ,WID,ISGRAV,IDGM,0,0,0)
C          Dalitz decay
          IF (MZIZ.GT.(AMGVSS+2*AME)) THEN
            WIDEE=WID*2*ALFAEM/(3*PI)*LOG(MZIZ/AME)
            CALL SSSAVE(ISZIZ,WIDEE,ISGRAV,IDE,-IDE,0,0)
          END IF
        END IF
        IF (MZIZ.GT.(AMZ+AMGVSS)) THEN
          WID=(2*(ZMIXSS(4,IZ)*STHW-ZMIXSS(3,IZ)*CTHW)**2+
     ,        (ZMIXSS(1,IZ)*SINBE-ZMIXSS(2,IZ)*COSBE)**2)*
     ,        (MZIZ**2-AMZ**2)**4/96./PI/MZIZ**3/(AMGVSS*AMPL)**2
          CALL SSSAVE(ISZIZ,WID,ISGRAV,IDZ,0,0,0)
        END IF
        IF (MZIZ.GT.(AMHL+AMGVSS)) THEN
          WID=(ZMIXSS(1,IZ)*COSA+ZMIXSS(2,IZ)*SINA)**2/6./
     ,        (AMGVSS*AMPL)**2*(MZIZ**2-AMHL**2)**4/16./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,WID,ISGRAV,ISHL,0,0,0)
        END IF
        IF (MZIZ.GT.(AMHH+AMGVSS)) THEN
          WID=(-ZMIXSS(1,IZ)*SINA+ZMIXSS(2,IZ)*COSA)**2/6./
     ,        (AMGVSS*AMPL)**2*(MZIZ**2-AMHH**2)**4/16./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,WID,ISGRAV,ISHH,0,0,0)
        END IF
        IF (MZIZ.GT.(AMHA+AMGVSS)) THEN
          WID=(ZMIXSS(1,IZ)*COSBE+ZMIXSS(2,IZ)*SINBE)**2/6./
     ,        (AMGVSS*AMPL)**2*(MZIZ**2-AMHA**2)**4/16./PI/MZIZ**3
          CALL SSSAVE(ISZIZ,WID,ISGRAV,ISHA,0,0,0)
        END IF
250   CONTINUE
C
C          Normalize zi branching ratios
C
      CALL SSNORM(ISZ1)
      CALL SSNORM(ISZ2)
      CALL SSNORM(ISZ3)
      CALL SSNORM(ISZ4)
C
      RETURN
      END
