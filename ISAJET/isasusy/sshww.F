#include "isajet/pilot.h"
      SUBROUTINE SSHWW
C-----------------------------------------------------------------------
C     Calculate HL, HH -> WW, ZZ, using either the on-shell matrix 
C     element if kinematically allowed or the WW* or ZZ* matrix
C     element from Eqn.(6) for Keung and Marciano (PRD. 84: 248).
C     For the latter, save the mode as W(Z) f fbar, and require that
C     MH > MW + 2 * MB.
C
C     Bisset's GBDCY
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstmp.inc"
#include "isajet/sstype.inc"
C
      EXTERNAL SSHWW1,SSHWW2
      DOUBLE PRECISION SSHWW1,SSHWW2
      DOUBLE PRECISION PI,SR2,G2,BETA,ALPHA,SW2,CW2,CAB2,SAB2,MW,MZ
     $,MH,COUPL,LOWER,UPPER,FWW1,FWW2,FWW3,FWW,DWID,FZZ
      DOUBLE PRECISION SSDINT,SSDLAM
      REAL WID
      REAL BRZN,BRZL,BRZU,BRZD,BRWL,BRWQ
      INTEGER IDHHA,IH
C          Hard wired Z branching ratios
      DATA BRZN,BRZL,BRZU,BRZD/.06839,.03442,.11792,.15191/
      DATA BRWL,BRWQ/.11111,.33333/
C
C          Mass matrix parameters
C
      PI=4*ATAN(1.D0)
      SR2=SQRT(2.D0)
      G2=4*PI*ALFAEM/SN2THW
      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
      SW2=SN2THW
      CW2=1.-SN2THW
      CAB2=(DCOS(ALPHA+BETA))**2
      SAB2=1.0-CAB2
      MW=AMW
      MZ=AMZ
C
C          WW* and ZZ* decays
C
      DO 100 IH=1,2
        IF(IH.EQ.1) THEN
          MH=AMHL
          IDHHA=ISHL
          COUPL=SAB2
        ELSE
          MH=AMHH
          IDHHA=ISHH
          COUPL=CAB2
        ENDIF
C          H -> W + W* -> W + f + fbar
        TMP(1)=MH
        IF(MH.GT.MW+2*AMBT.AND.MH.LE.2*MW) THEN
          LOWER=2*MW/MH
          UPPER=1+MW**2/MH**2
          IF (LOWER.LT.0.998D0) THEN
            IF (UPPER.LE.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,0.998D0)
              FWW2=SSDINT(0.998D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2
            ELSEIF(UPPER.GT.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,0.998D0)
              FWW2=SSDINT(0.998D0,SSHWW1,1.001D0)
              FWW3=SSDINT(1.001D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2+FWW3
            ENDIF
          ELSE IF (0.998D0.LT.LOWER.AND.LOWER.LT.1.001D0) THEN
            IF (UPPER.LE.1.001D0) THEN
              FWW=SSDINT(LOWER,SSHWW1,UPPER)
            ELSEIF(UPPER.GT.1.001D0) THEN
              FWW1=SSDINT(LOWER,SSHWW1,1.001D0)
              FWW2=SSDINT(1.001D0,SSHWW1,UPPER)
              FWW=FWW1+FWW2
            ENDIF
          ELSE IF (LOWER.GT.1.001D0) THEN
            FWW=SSDINT(LOWER,SSHWW1,UPPER)
          END IF
          DWID=3*(G2**2)*MH*FWW/(512.0*PI**3)
          WID=DWID*COUPL
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDE,-IDNE,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDMU,-IDNM,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,IDW,IDTAU,-IDNT,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,IDW,-IDUP,IDDN,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,IDW,-IDCH,IDST,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDE,IDNE,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDMU,IDNM,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWL*WID,-IDW,-IDTAU,IDNT,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,-IDW,IDUP,-IDDN,0,0)
          CALL SSSAVE(IDHHA,0.5*BRWQ*WID,-IDW,IDCH,-IDST,0,0)
        ENDIF
C          H -> Z + Z* -> Z + f + fbar
        IF(MH.GT.MZ+2*AMBT.AND.MH.LE.2*MZ) THEN
          LOWER=2*MZ/MH
          UPPER=1+MZ**2/MH**2
          FZZ=SSDINT(LOWER,SSHWW2,UPPER)
          DWID=7.0-40*SW2/3+160*SW2**2/9
          DWID=DWID/CW2**2             
          DWID=DWID*G2**2*MH*FZZ/(2048*PI**3)
          WID=DWID*COUPL
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNE,-IDNE,0,0)
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNM,-IDNM,0,0)
          CALL SSSAVE(IDHHA,BRZN*WID,IDZ,IDNT,-IDNT,0,0)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDE,-IDE,0,0)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDMU,-IDMU,0,0)
          CALL SSSAVE(IDHHA,BRZL*WID,IDZ,IDTAU,-IDTAU,0,0)
          CALL SSSAVE(IDHHA,BRZU*WID,IDZ,IDUP,-IDUP,0,0)
          CALL SSSAVE(IDHHA,BRZU*WID,IDZ,IDCH,-IDCH,0,0)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDDN,-IDDN,0,0)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDST,-IDST,0,0)
          CALL SSSAVE(IDHHA,BRZD*WID,IDZ,IDBT,-IDBT,0,0)
        ENDIF
100   CONTINUE
C
C          HH -> WW, ZZ
C          If these are allowed, the WW* and ZZ* are not.
C
      MH=AMHH
      IF(MH.GT.2*MW) THEN
        DWID=3+(MH/MW)**4/4-(MH/MW)**2
        DWID=DWID*G2*CAB2*MW**2/(16.0*PI*MH**3)
        WID=DWID*SQRT(SSDLAM(MH**2,MW**2,MW**2))
        CALL SSSAVE(ISHH,WID,IDW,-IDW,0,0,0)
      ENDIF
      IF(MH.GT.2*MZ) THEN
        DWID=3+(MH/MZ)**4/4-(MH/MZ)**2
        DWID=DWID*G2*CAB2*MW**2/(16.0*PI*MH**3)/(2.0*CW2**2)
        WID=DWID*SQRT(SSDLAM(MH**2,MZ**2,MZ**2))
        CALL SSSAVE(ISHH,WID,IDZ,IDZ,0,0,0)
      ENDIF
C
      RETURN
      END
