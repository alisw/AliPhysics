#include "isajet/pilot.h"
      LOGICAL FUNCTION DECTAU(IP,NADD,MEIP,IDABS,PREST)
C
C          Compute matrix elements for polarized tau decay.
C          Polarization determined by tau parent.
C          Auxiliary routine for DECAY. 
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/wcon.inc"
#include "isajet/partcl.inc"
#include "isajet/dkytab.inc"
#include "isajet/const.inc"
#include "isajet/pjets.inc"
#include "isajet/keys.inc"
#include "isajet/xmssm.inc"
#include "isajet/sspols.inc"
#include "isajet/primar.inc"
C
      REAL PREST(4,6),WT,TAUHEL,S12,S12MAX,PIP,CTHNU,PSUM(4),AMV2,WT1
      REAL DOT,DOT3,RANF,Z
      INTEGER IP,NADD,IDABS(5),IPAR,IDPAR,JET,INU,I,K,I1,I2,IDSIB
      INTEGER IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX,IDIP
      INTEGER MEIP,IPX,IP1,IP2
C
      DOT(I1,I2)=PREST(4,I1)*PREST(4,I2)-PREST(1,I1)*PREST(1,I2)
     $-PREST(2,I1)*PREST(2,I2)-PREST(3,I1)*PREST(3,I2)
      DOT3(I1,I2)=PREST(1,I1)*PREST(1,I2)+PREST(2,I1)*PREST(2,I2)
     $+PREST(3,I1)*PREST(3,I2)
C
      IDIP=IDENT(IP)
      DECTAU=.TRUE.
      IF(IABS(IDIP).NE.16) GO TO 999
C
C          Use PREST(K,6) for spin vector
C
      PIP=SQRT(PPTCL(1,IP)**2+PPTCL(2,IP)**2+PPTCL(3,IP)**2)
      DO 100 K=1,3
        PREST(K,6)=PPTCL(K,IP)/PIP
100   CONTINUE
      PREST(4,6)=0.
C
C          Take helicity TAUHEL=0 unless TAU parent is TP, W+-, H+-,
C          or some SUSY particles.
C          Take account of 1-particle decays!
C

      IPX=IP
      TAUHEL=0.
      IPAR=0
      IDPAR=0
110   IF(IORIG(IPX).GT.0) THEN
        IPAR=MOD(IORIG(IPX),IPACK)
        IDPAR=IDENT(IPAR)
        IF(IDPAR.EQ.IDIP) THEN
          IP1=IDCAY(IPAR)/IPACK
          IP2=MOD(IDCAY(IPAR),IPACK)
          IF(IP1.EQ.IP2) THEN
            IPX=IPAR
            GO TO 110
          ENDIF
        ENDIF
        IDPAR=IABS(IDPAR)
        IDSIB=0
C          W/top parent
        IF(IDPAR.GT.100.AND.MOD(IDPAR/10,10).GE.6) THEN
          TAUHEL=-1.
        ELSEIF(IDPAR.EQ.80) THEN
          TAUHEL=-1.
C          Charged Higgs parent
        ELSEIF(IDPAR.EQ.86) THEN
          TAUHEL=+1.
C          SUSY parent - polarization also depends on sibling IDSIB
        ELSEIF(GOMSSM.AND.IDPAR.GT.20.AND.IDPAR.LT.80) THEN
          I1=IDCAY(IPAR)/IPACK
          I2=MOD(IDCAY(IPAR),IPACK)
          DO 120 I=I1,I2
            IF(IABS(IDENT(I)).GT.20.AND.IABS(IDENT(I)).LT.80)
     $      IDSIB=IABS(IDENT(I))
120       CONTINUE
          IF (IDPAR.EQ.35) THEN
            TAUHEL=-1.
          ELSEIF (IDPAR.EQ.36) THEN
            IF (IDSIB.EQ.30) TAUHEL=PTAU1(1)
            IF (IDSIB.EQ.40) TAUHEL=PTAU1(2)
            IF (IDSIB.EQ.50) TAUHEL=PTAU1(3)
            IF (IDSIB.EQ.60) TAUHEL=PTAU1(4)
          ELSEIF (IDPAR.EQ.56) THEN
            IF (IDSIB.EQ.30) TAUHEL=PTAU2(1)
            IF (IDSIB.EQ.40) TAUHEL=PTAU2(2)
            IF (IDSIB.EQ.50) TAUHEL=PTAU2(3)
            IF (IDSIB.EQ.60) TAUHEL=PTAU2(4)
          ELSEIF (IDPAR.EQ.39) THEN
            IF(IDSIB.EQ.35) TAUHEL=-1.
            IF(IDSIB.EQ.30) TAUHEL=PTAUWZ
          ELSEIF (IDPAR.EQ.49.AND.IDSIB.EQ.35) THEN
            TAUHEL=-1.
          ELSEIF (IDPAR.EQ.40) THEN
            IF(IDSIB.EQ.36) TAUHEL=PTAUZ2(1)
            IF(IDSIB.EQ.56) TAUHEL=PTAUZ2(2)
            IF(IDSIB.EQ.30) TAUHEL=PTAUZZ
          ELSEIF (IDPAR.EQ.50) THEN
            IF(IDSIB.EQ.36) TAUHEL=PTAUZ3(1)
            IF(IDSIB.EQ.56) TAUHEL=PTAUZ3(2)
          ELSEIF (IDPAR.EQ.60) THEN 
            IF(IDSIB.EQ.36) TAUHEL=PTAUZ4(1)
            IF(IDSIB.EQ.56) TAUHEL=PTAUZ4(2)
          ENDIF
        END IF
      ELSE
        IF(KEYS(3)) THEN
          IF(IABS(IDENTW).EQ.80) TAUHEL=-1.
        ELSE
          JET=IABS(IORIG(IP))/IPACK
          IF(JET.GT.0.AND.JET.LE.NJET) THEN
            IF(IDJETS(JET).EQ.80) TAUHEL=-1.
          ENDIF
        ENDIF
      ENDIF
C
C          Leptonic decays. DECTAU is always called for TAU- decay
C          products, so selection is independent of IDENT(IP).
C
      IF(MEIP.EQ.5) THEN
        IF(IDENT(NPTCL+1).LT.0) THEN
          WT=PPTCL(5,IP)*(PREST(4,1)-TAUHEL*DOT(1,6))*DOT(2,3)
        ELSEIF(IDENT(NPTCL+2).LT.0) THEN
          WT=PPTCL(5,IP)*(PREST(4,2)-TAUHEL*DOT(2,6))*DOT(1,3)
        ELSE
          WT=PPTCL(5,IP)*(PREST(4,3)-TAUHEL*DOT(3,6))*DOT(1,2)
        ENDIF
        IF(WT.LT.RANF()*PPTCL(5,IP)**4/8.) THEN
          DECTAU=.FALSE.
        ELSE
          DECTAU=.TRUE.
        ENDIF
        RETURN
C
C          Decay to PI + NUT, K + NUT
C
      ELSEIF(MEIP.EQ.6) THEN
        INU=1
        IF(IDABS(2).EQ.15) INU=2
        CTHNU=DOT3(INU,6)/SQRT(DOT3(INU,INU))
        WT=1.-TAUHEL*CTHNU
        IF(WT.LT.RANF()*2.) THEN
          DECTAU=.FALSE.
        ELSE
          DECTAU=.TRUE.
        ENDIF
        RETURN
C
C          Decay to RHO + NUT, A1 + NUT, K* + NUT
C
      ELSEIF(MEIP.EQ.7) THEN
        DO 210 I=1,NADD
210     IF(IDABS(I).EQ.15) INU=I
        DO 220 K=1,4
          PSUM(K)=0.
          DO 221 I=1,NADD
            IF(I.EQ.INU) GO TO 221
            PSUM(K)=PSUM(K)+PREST(K,I)
221       CONTINUE
220     CONTINUE
        AMV2=PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2
        WT1=2.*AMV2/(2.*AMV2+PPTCL(5,IP)**2)
        CTHNU=DOT3(INU,6)/SQRT(DOT3(INU,INU))
        WT=WT1*(1.+TAUHEL*CTHNU)+(1.-WT1)*(1-TAUHEL*CTHNU)
        IF(WT.LT.RANF()*2.) THEN
          DECTAU=.FALSE.
        ELSE
          DECTAU=.TRUE.
        ENDIF
        RETURN
C
C          Ignore matrix element for all other decays
C
      ELSE
        DECTAU=.TRUE.
        RETURN
      ENDIF
C
C          Error
C
999   CALL PRTEVT(0)
      WRITE(ITLIS,99999) IP
99999 FORMAT(//5X,'ERROR IN DECTAU FOR PARTICLE',I6)
      END
