#include "isajet/pilot.h"
C---------------------------------------------------------------
      SUBROUTINE SUGMAS(G0,ILOOP,IMODEL)
C---------------------------------------------------------------
C
C     Compute tree level sparticle masses; output to MSS, XISAIN
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sugpas.inc"
#include "isajet/sugxin.inc"
#include "isajet/sugmg.inc"
      REAL MSB1,MSB2,MST1,MST2
      REAL G0(29)
      REAL SUGMFN,SUALFS,SSPOLE,MHP,MGLMGL,MHPS,
     $RDEL,ASMGL,DELHPS,M1S,M2S,FNB,FCN,
     $MB,FNT,MT,MW,TANB,BETA,COSB,COTB,SINB,MZ,COS2B,
     $PI,T2S,G,ATAU,MSSS,AT,AB,BRKT,B2S,T1S,TERM,B1S,Q,
     $MBQ,MTAMZ,MTQ,FNL,MSL1,MSL2,ASMB,MBMB,ASMT,MTMT
      REAL AA,BB,CC,DA,DB,DC,L1,L2,EVAL1,RL1,RL2
      DOUBLE PRECISION SSMQCD
      INTEGER IALLOW,ILOOP,MHLNEG,MHCNEG,IMODEL
C
C          Statement function
C
      SUGMFN(Q)=Q**2*(LOG(Q**2/HIGFRZ**2)-1.)
C
      PI=4.*ATAN(1.)
      XW=.232
      G=G2
      TANB=XTANB
      MT=AMT
      MZ=AMZ
      MW=AMW
      AMTP=MT
      BETA=ATAN(TANB)
      COTB=1./TANB
      SINB=SIN(BETA)
      COSB=COS(BETA)
      SIN2B=SIN(2*BETA)
      COS2B=COS(2*BETA)
      AT=G0(12)
      AB=G0(11)
      ATAU=G0(10)
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(HIGFRZ))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(HIGFRZ))
      MTAMZ=FTAMZ*COSB*VEV
C
C          Compute some masses from RGE solution to prepare for SSMASS,
C          which computes the rest.
C
      MSSS=G0(19)+AMUP**2+(.5-2*XW/3.)*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
C          Squark and slepton masses
      MSS(2)=SQRT(MSSS)
      MSS(3)=SQRT(G0(18)+AMUP**2+2./3.*XW*MZ**2*COS2B)
      MSS(4)=SQRT(G0(19)+AMDN**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(5)=SQRT(G0(17)+AMDN**2-1./3.*XW*MZ**2*COS2B)
      MSS(6)=SQRT(G0(19)+AMST**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(7)=SQRT(G0(17)+AMST**2-1./3.*XW*MZ**2*COS2B)
      MSS(8)=SQRT(G0(19)+AMCH**2+(.5-2*XW/3.)*MZ**2*COS2B)
      MSS(9)=SQRT(G0(18)+AMCH**2+2./3.*XW*MZ**2*COS2B)
      BRKT=(.5*(G0(24)-G0(22))-COS2B*(4*MW**2-MZ**2)/12.)**2+
     $       MBQ**2*(AB-MU*TANB)**2
      TERM=.5*(G0(24)+G0(22))+MBQ**2-MZ**2*COS2B/4.
      B1S=TERM-SQRT(BRKT)
      B2S=TERM+SQRT(BRKT)
      MSS(10)=SQRT(MAX(0.,B1S))
      MSS(11)=SQRT(MAX(0.,B2S))
      BRKT=(.5*(G0(24)-G0(23))+COS2B*(8*MW**2-5*MZ**2)/12.)**2+
     $       MTQ**2*(AT-MU*COTB)**2
      TERM=.5*(G0(24)+G0(23))+MTQ**2+MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0..OR.B1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(12)=SQRT(MAX(0.,T1S))
      MSS(13)=SQRT(MAX(0.,T2S))
      MSSS=G0(16)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(14)=SQRT(MSSS)
      MSS(15)=MSS(14)
      MSSS=G0(21)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(16)=SQRT(MSSS)
      MSS(17)=SQRT(G0(16)+AME**2-.5*(2*MW**2-MZ**2)*COS2B)
      MSS(18)=SQRT(G0(15)+AME**2+(MW**2-MZ**2)*COS2B)
      MSS(19)=SQRT(G0(16)+AMMU**2-.5*(2*MW**2-MZ**2)*COS2B)
      MSS(20)=SQRT(G0(15)+AMMU**2+(MW**2-MZ**2)*COS2B)
      BRKT=(.5*(G0(21)-G0(20))-COS2B*(4*MW**2-3*MZ**2)/4.)**2+
     $       MTAMZ**2*(ATAU-MU*TANB)**2
      TERM=.5*(G0(21)+G0(20))+MTAMZ**2-MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(21)=SQRT(MAX(0.,T1S))
      MSS(22)=SQRT(MAX(0.,T2S))
C          A0 mass
      M1S=MU**2+G0(13)
      M2S=MU**2+G0(14)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MST1=MSS(12)
      MST2=MSS(13)
      MSL1=MSS(21)
      MSL2=MSS(22)
      MB=AMBT
      FNT=(SUGMFN(MST2)-SUGMFN(MST1))/(MST2**2-MST1**2)
     $*AT*MTQ**2/SINB**2
      FNB=(SUGMFN(MSB2)-SUGMFN(MSB1))/(MSB2**2-MSB1**2)
     $*AB*MBQ**2/COSB**2
      FNL=(SUGMFN(MSL2)-SUGMFN(MSL1))/(MSL2**2-MSL1**2)
     $*ATAU*MTAMZ**2/COSB**2
      FCN=FNT+FNB+FNL/3.
      DELHPS=3*G0(2)**2*MU*(COTB+TANB)/32./PI**2/MW**2*FCN
      RDEL=SQRT(ABS(DELHPS))
C          Tree level mhp not needed at this point so fix if negative
      IF (ILOOP.EQ.0) THEN
        MHPS=M1S+M2S
        IF (MHPS.LT.0.) MHPS=0.
      ELSE
        MHPS=B*MU*(COTB+TANB)+DELHPS
        IF (MHPS.LT.0.) THEN
          NOGOOD=3
          MHPS=AMZ**2
        END IF
      END IF
      MHP=SQRT(MHPS)
      MSS(31)=MHP
C     APPLY XERXES' TEST FOR PROPER POTENTIAL SHAPE AT THE ORIGIN
C     REMOVE THIS CONSTRAINT ON 4/7/00
      IF (ILOOP.EQ.1) THEN
      L1=MIN(G0(24),G0(23))
      L2=MAX(G0(24),G0(23))
      RL1=SQRT(L1)
      RL2=SQRT(L2)
      DA=3*G0(6)**2*AT**2/ABS(G0(24)-G0(23))/16./PI**2*
     $(-SUGMFN(RL1)+SUGMFN(RL2))
      DB=3*G0(6)**2/16./PI**2*
     $(SUGMFN(RL1)*(1.-AT**2/ABS(G0(24)-G0(23)))+SUGMFN(RL2)*
     $(1.+AT**2/ABS(G0(24)-G0(23))))
      DC=-3*G0(6)**2*AT*MU/ABS(G0(24)-G0(23))/16./PI**2*
     $(-SUGMFN(RL1)+SUGMFN(RL2))
      AA=M1S+DA
      BB=M2S+DB
      CC=-B*MU+DC
      EVAL1=((AA+BB)-SQRT((AA+BB)**2-4*(AA*BB-CC*CC)))/2.
C      IF (EVAL1.GE.0) THEN
C        NOGOOD=7
C      END IF
      END IF
C
C          Initialize SUSY parameters in /SSPAR/:
C
      AMGLSS=G0(9)
      AMULSS=MSS(2)
      AMURSS=MSS(3)
      AMDLSS=MSS(4)
      AMDRSS=MSS(5)
      AMSLSS=MSS(6)
      AMSRSS=MSS(7)
      AMCLSS=MSS(8)
      AMCRSS=MSS(9)
      AMN1SS=MSS(16)
      AMN2SS=MSS(16)
      AMN3SS=MSS(16)
      AMELSS=MSS(17)
      AMERSS=MSS(18)
      AMMLSS=MSS(19)
      AMMRSS=MSS(20)
      TWOM1=-MU
      RV2V1=1./TANB
      AMTLSS=SQRT(G0(24))
      AMTRSS=SQRT(G0(23))
      AMBLSS=SQRT(G0(24))
      AMBRSS=SQRT(G0(22))
      AMLLSS=SQRT(G0(21))
      AMLRSS=SQRT(G0(20))
      AAT=G0(12)
      AAB=G0(11)
      AAL=G0(10)
      AMHA=MHP
C
C          Use SSMASS to diagonalize neutralino and chargino mass
C          matrices and calculate Higgs masses.
C
      MHLNEG=0
      MHCNEG=0
      CALL SSMASS(G0(7),G0(8),IALLOW,ILOOP,MHLNEG,MHCNEG,IMODEL)
      IF(MHLNEG.EQ.1.OR.MHCNEG.EQ.1) THEN
        NOGOOD=8
      ENDIF
      IF(IALLOW.NE.0) THEN
        NOGOOD=5
        GO TO 100
      ENDIF
C
C          Save results also in MSS
C
      MSS(23)=AMZ1SS
      MSS(24)=AMZ2SS
      MSS(25)=AMZ3SS
      MSS(26)=AMZ4SS
      MSS(27)=AMW1SS
      MSS(28)=AMW2SS
      MSS(29)=AMHL
      MSS(30)=AMHH
      MSS(31)=AMHA
      MSS(32)=AMHC
C          Gluino pole mass
      MGLMGL=G0(9)
      ASMGL=SUALFS(MGLMGL**2,.36,MT,3)
      MSS(1)=SSPOLE(MGLMGL,MGLMGL**2,ASMGL)
      AMGLSS=MSS(1)
C
100   RETURN
      END
