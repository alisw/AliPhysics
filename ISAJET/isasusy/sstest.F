#include "isajet/pilot.h"
      SUBROUTINE SSTEST(IALLOW)
C
C          Test MSSM parameters against existing bounds on SUSY from
C          LEP and SLC:
C          IALLOW = 1    Z1 is not LSP
C          IALLOW = 2    Gamma(Z -> Z1SS Z1SS) < GAMINV
C          IALLOW = 4    Z -> charginos allowed
C          IALLOW = 8    BF(Z -> Z1SS Z2SS)>10^5
C          IALLOW = 16   Z -> squarks, sleptons
C          IALLOW = 32   BR(Z -> Z* HL0) < B(Z -> Z* H(M=MHSM))
C          IALLOW = 64   BR(Z -> HL0 HA0) > 0
C          IALLOW = 128  M(H+) > M(Z)/2
C          where GAMINV is the present bound on the invisible width,
C          and MHSM is the lower bound on the standard Higgs mass.
C
C          Bounds on the other modes are only approximate, but the
C          error in the allowed region of masses must be tiny. 
C          Updated by H. Baer on 5/25/95
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
#include "isajet/sstype.inc"
#include "isajet/sugmg.inc"
#include "isajet/xmssm.inc"
C
      INTEGER IALLOW
      EXTERNAL SSZHX
      REAL MHSM,GAMINV,PI,SR2,G,GP,MZ,MZ1,MZ2,MZ3,MZ4,MW1,MW2,
     $TANB,BETA,COS2B,SIN2B,VS,V,VP,FT,MHL,ALPHA,SUSYCC,
     $GAMSS,W11,GZ1Z1,GAMSM,SSXINT,SSXLAM,COS2W,
     $W12,GZ1Z2,DGAMZ,BFZ,BFZ1Z2
C
C          Current bounds
      DATA MHSM/64./,GAMINV/.0043/,DGAMZ/.0115/,BFZ/1.E-5/
C
C          Initialize
C
      IALLOW=0
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      COS2W=SQRT(1.-SN2THW)
      MZ=AMZ
      MZ1=ABS(AMZ1SS)
      MZ2=ABS(AMZ2SS)
      MZ3=ABS(AMZ3SS)
      MZ4=ABS(AMZ4SS)
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
C
C          Check that Z1SS is LSP
C
      IF(MZ1.GT.MW1.OR.MZ1.GT.AMGLSS.OR.MZ1.GT.AMULSS
     $.OR.MZ1.GT.AMERSS.OR.MZ1.GT.AMELSS.OR.MZ1.GT.AMN1SS
     $.OR.MZ1.GT.AMB1SS.OR.MZ1.GT.AMT1SS.OR.MZ1.GT.AML1SS) THEN
        IALLOW=IALLOW+1
      ENDIF
C
C          Z -> Z1SS + Z1SS
C
      IF (AMZ.GT.2*MZ1) THEN
        W11=SQRT(G*G+GP*GP)
     $  *(ZMIXSS(1,1)*ZMIXSS(1,1)-ZMIXSS(2,1)*ZMIXSS(2,1))/4.
        GZ1Z1=SQRT(SSXLAM(MZ**2,MZ1**2,MZ1**2))/12./PI/MZ**3*W11**2
     $  *(2*MZ**2-MZ1**2-MZ1**2-(MZ1**2-MZ1**2)**2/MZ**2
     $  -6*MZ1*MZ1*SIGN(1.,AMZ1SS*AMZ1SS))
        IF(GZ1Z1.GT.GAMINV) THEN
          IALLOW=IALLOW+2
        ENDIF
      ENDIF
C
C          Check for other allowed visible modes modes
C
      IF(AMZ.GT.2*MW1) THEN
        IALLOW=IALLOW+4
      ENDIF
C
C          Check funny Z branching fractions
C
      BFZ1Z2=0.
      IF (AMZ.GT.MZ1+MZ2) THEN
        W12=SQRT(G*G+GP*GP)
     $  *(ZMIXSS(1,1)*ZMIXSS(1,2)-ZMIXSS(2,1)*ZMIXSS(2,2))/4.
        GZ1Z2=SQRT(SSXLAM(MZ**2,MZ1**2,MZ2**2))/6./PI/MZ**3*W12**2
     $  *(2*MZ**2-MZ1**2-MZ2**2-(MZ1**2-MZ2**2)**2/MZ**2
     $  -6*MZ1*MZ2*SIGN(1.,AMZ1SS*AMZ2SS))
        BFZ1Z2=GZ1Z2/GAMZ
      END IF
      IF(BFZ1Z2.GT.BFZ) THEN
        IALLOW=IALLOW+8
      ENDIF
C
      IF(AMZ.GT.2*AMULSS.OR.AMZ.GT.2*AMELSS.OR.AMZ.GT.2*AMERSS
     $.OR.AMZ.GT.2*AMN1SS.OR.AMZ.GT.2*AMB1SS.OR.AMZ.GT.2*AMT1SS)THEN
        IALLOW=IALLOW+16
      ENDIF
C
C          Z -> Higgs modes
C
      TMP(1)=MHSM
      GAMSM=SSXINT(2*MHSM/MZ,SSZHX,(1.+MHSM**2/MZ**2))  
C          Z -> hl0 x
      IF(AMZ.GT.AMHL) THEN
        TANB=1./RV2V1
        BETA=ATAN(TANB)
        COS2B=COS(2*BETA)
        SIN2B=SIN(2*BETA)
        VS=2*AMW**2/G**2/(1.+RV2V1**2)
        V=SQRT(VS)
        VP=RV2V1*V
        FT=G*AMTP/SR2/AMW/V*SQRT(V**2+VP**2)
        MHL=AMHL
        ALPHA=ALFAH
        SUSYCC=SIN(ALPHA+BETA)
        TMP(1)=MHL
        GAMSS=SSXINT(2*MHL/AMZ,SSZHX,(1.+MHL**2/AMZ**2))*SUSYCC**2
        IF(GAMSS.GE.GAMSM) IALLOW=IALLOW+32
      ENDIF
C          Z -> hl0 ha0
      IF (AMZ.GT.(AMHL+AMHA)) THEN
        IALLOW=IALLOW+64
      ENDIF
C          Z -> H+ H-
      IF(AMZ.GT.2*AMHC) THEN
        IALLOW=IALLOW+128
      ENDIF
C
      RETURN
      END
