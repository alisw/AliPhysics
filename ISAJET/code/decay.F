#include "isajet/pilot.h"
      SUBROUTINE DECAY(IP)
C
C          Decay particle IP from /PARTCL/ using /DKYTAB/ branching
C          ratios and add decay products to /PARTCL/ with IORIG=IP.
C          Forced decay modes are flagged by LOOK<0.
C
C          Auxiliary routines:
C          DECPS1: generate masses for phase space
C          DECPS2: generate 2-body decays and boosts for phase space
C          DECVA:  V-A matrix elements
C          DECTAU: tau decay matrix elements with polarization
C          DECSS3: 3-body SUSY matrix element using /DKYSS3/
C          DECJET: Hadronize partons from decay.
C
C          Matrix element for Dalitz decays and W mass for TP -> W BT
C          are generated explicitly. W width is included.
C
C          Requirements for decay modes:
C          (1) For Dalitz decays, particle 1 must be GM.
C          (2) For V-A quark or lepton decays, particles 1 and 2 must
C              be from (virtual) W.
C          (3) For any decay into quarks, they must appear last.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/wcon.inc"
#include "isajet/partcl.inc"
#include "isajet/dkytab.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/const.inc"
#include "isajet/primar.inc"
#include "isajet/idrun.inc"
#include "isajet/force.inc"
#include "isajet/sstype.inc"
#include "isajet/dkyss3.inc"
C
      REAL PGEN(5,5),BETA(3),REDUCE(5),WPROP,Z,TRY,RANF,AMASS,TWOME
      REAL PSUM(5),SUM,PREST(4,6),DOT,PCM
      REAL AMEE,REE,WTEE,SWAP,WT,A,B,C,GAMMA
      REAL SMAX,SMIN,SVAL,TANMAX,TANMIN,TANVAL
      LOGICAL WDECAY,DECVA,DECTAU,DECJET
      INTEGER IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX,IPOINT,ID1,I1,I2
      INTEGER NADD,NSTART,NEW,NADD1,J,IP,I,IDABS(5)
      INTEGER K,JETIP,IDANTI,NPASS,MEIP,MEA
      REAL DBLPCM,DECSS3,VAL
C
      DATA REDUCE/1.,1.,2.,5.,15./
      DATA PSUM/5*0./
      DATA TWOME/1.022006E-3/
      DATA PREST/24*0./
C
C          Function definitions.
C          Use double precision for PCM on 32-bit machines
C
#if defined(CERNLIB_SINGLE)
      PCM(A,B,C)=SQRT((A**2-B**2-C**2)**2-(2.*B*C)**2)/(2.*A)
#endif
#if defined(CERNLIB_DOUBLE)
      PCM(A,B,C)=DBLPCM(A,B,C)
#endif
      DOT(I1,I2)=PREST(4,I1)*PREST(4,I2)-PREST(1,I1)*PREST(1,I2)
     $-PREST(2,I1)*PREST(2,I2)-PREST(3,I1)*PREST(3,I2)
C          Charged W propagator.
      WPROP(Z)=(Z-WMASS(2)**2)**2+(WMASS(2)*WGAM(2))**2
C----------------------------------------------------------------------
C          Select decay mode. Note IDENT(NPTCL+1)...IDENT(NPTCL+5)
C          are always defined even if zero.
C----------------------------------------------------------------------
      IF(IDCAY(IP).NE.0) RETURN
      IDLV1=IDENT(IP)
      CALL FLAVOR(IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX)
C          FLAVOR returns 0 for quark, but want IFL3=6 for top
      IF(IABS(IDLV1).LT.10) IFL3=IDLV1
      NPASS=0
1     CONTINUE
      NPASS=NPASS+1
      WDECAY=.FALSE.
      IF(NPASS.GT.NTRIES) GO TO 9998
      IPOINT=LOOK(INDEX)
      IF(IPOINT.EQ.0) RETURN
C          IPOINT<0 flags a forced decay.
      IF(IPOINT.LT.0) THEN
        I=1
        IF(IDENT(IP).LT.0) I=2
        IPOINT=LOOK2(I,IABS(IPOINT))
      ENDIF
C
C          Select mode.
C
      TRY=RANF()
      IPOINT=IPOINT-1
100   IPOINT=IPOINT+1
      IF(TRY.GT.CBR(IPOINT)) GO TO 100
      NADD=0
      SUM=0.
      NSTART=NPTCL+1
      IF(NPTCL+5.GT.MXPTCL) GO TO 9999
C
C          Set up masses and IDENT codes.
C
      MEIP=MELEM(IPOINT)
      DO 110 I=1,5
        NEW=NPTCL+I
        IDENT(NEW)=MODE(I,IPOINT)
        IDABS(I)=IABS(IDENT(NEW))
        IF(MODE(I,IPOINT).EQ.0) GO TO 110
        NADD=NADD+1
        IDLV1=IDENT(NEW)
        PPTCL(5,NEW)=AMASS(IDLV1)
        SUM=SUM+PPTCL(5,NEW)
110   CONTINUE
      NADD1=NADD-1
      DO 120 J=1,5
        PGEN(J,1)=PPTCL(J,IP)
120   CONTINUE
      PGEN(5,NADD)=PPTCL(5,NPTCL+NADD)
C----------------------------------------------------------------------
C          Carry out appropriate decay
C----------------------------------------------------------------------
C
C          1-body decays.
C
      IF(NADD.EQ.1) THEN
        DO 200 J=1,5
          PPTCL(J,NPTCL+1)=PPTCL(J,IP)
200     CONTINUE
        GO TO 300
      ENDIF
C
C          2-body phase space decays
C
      IF(NADD.EQ.2.AND.MEIP.EQ.0) THEN
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          N-body phase space decays
C
      IF(NADD.GT.2.AND.MEIP.EQ.0) THEN
        CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          Dalitz decays
C
      IF(NADD.EQ.3.AND.MEIP.EQ.1) THEN
210     AMEE=TWOME*(PPTCL(5,IP)/TWOME)**RANF()
        REE=(TWOME/AMEE)**2
        WTEE=(1.-(AMEE/PPTCL(5,IP))**2)**3*SQRT(1.-REE)*(1.+.5*REE)
        IF(WTEE.LT.RANF()) GO TO 210
        PGEN(5,2)=AMEE
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        GO TO 300
      ENDIF
C
C          omega/phi decays (for reasons lost in history...)
C
      IF(NADD.EQ.3.AND.MEIP.EQ.2) THEN
220     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        WT=(PPTCL(5,NPTCL+1)*PPTCL(5,NPTCL+2)*PPTCL(5,NPTCL+3))**2
     $  -(PPTCL(5,NPTCL+1)*DOT(2,3))**2
     $  -(PPTCL(5,NPTCL+2)*DOT(1,3))**2
     $  -(PPTCL(5,NPTCL+3)*DOT(1,2))**2
     $  +2.*DOT(1,2)*DOT(2,3)*DOT(1,3)
        IF(WT.LT.RANF()*PPTCL(5,IP)**6/108.) GO TO 220
        GO TO 300
      ENDIF
C
C          V-A decays
C
      IF(NADD.EQ.3.AND.MEIP.EQ.3) THEN
230     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECVA(IP,NADD,IDABS,PREST)) GO TO 230
        GO TO 300
      ENDIF
C
C          Top decays
C          Generate mass for TP -> W BT with Breit-Wigner. 
C          W couples to 1+2 so swap 1<->3. Then m2+m3 < m < m0-m1.
C
      IF(NADD.EQ.3.AND.MEIP.EQ.4) THEN
        WDECAY=.TRUE.
        SWAP=PPTCL(5,NPTCL+1)
        PPTCL(5,NPTCL+1)=PPTCL(5,NPTCL+3)
        PPTCL(5,NPTCL+3)=SWAP
        SMAX=(PPTCL(5,IP)-PPTCL(5,NPTCL+1))**2
        SMIN=(PPTCL(5,NPTCL+2)+PPTCL(5,NPTCL+3))**2
        TANMAX=ATAN((SMAX-WMASS(2)**2)/(WMASS(2)*WGAM(2)))
        TANMIN=ATAN((SMIN-WMASS(2)**2)/(WMASS(2)*WGAM(2)))
240     TANVAL=RANF()*(TANMAX-TANMIN)+TANMIN
        SVAL=WMASS(2)**2+WMASS(2)*WGAM(2)*TAN(TANVAL)
        IF(SVAL.LT.SMIN.OR.SVAL.GT.SMAX) GO TO 240
        PGEN(5,2)=SQRT(SVAL)
        PGEN(5,3)=PPTCL(5,NPTCL+3)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECVA(IP,NADD,IDABS,PREST)) GO TO 240
        DO 241 K=1,5
          SWAP=PPTCL(K,NPTCL+1)
          PPTCL(K,NPTCL+1)=PPTCL(K,NPTCL+3)
          PPTCL(K,NPTCL+3)=SWAP
241     CONTINUE
        PGEN(5,3)=PPTCL(5,NPTCL+3)
        DO 242 K=1,4
          SWAP=PREST(K,1)
          PREST(K,1)=PREST(K,3)
          PREST(K,3)=SWAP
242     CONTINUE
        GO TO 300
      ENDIF
C
C          TAU decays. These are special because they take polarization
C          into account.
C
      IF(MEIP.EQ.5.OR.MEIP.EQ.6.OR.MEIP.EQ.7) THEN
250     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        IF(.NOT.DECTAU(IP,NADD,MEIP,IDABS,PREST)) GO TO 250
        GO TO 300
      ENDIF
C
C          3-body SUSY decays
C
      IF(MEIP.LT.0.AND.NADD.EQ.3) THEN
        MEA=IABS(MEIP)
        IF(WTSS3(MEA).LE.0) THEN
          DO 260 I=1,1000
            CALL DECPS1(IP,NADD,PGEN)
            CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
            VAL=DECSS3(IP,MEA)
            WTSS3(MEA)=MAX(WTSS3(MEA),VAL)
260       CONTINUE
          IF(WTSS3(MEA).LE.0) GO TO 9998
        ENDIF
261     CALL DECPS1(IP,NADD,PGEN)
        CALL DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
        VAL=DECSS3(IP,MEA)
        WTSS3(MEA)=MAX(WTSS3(MEA),VAL)
        IF(VAL.LT.WTSS3(MEA)*RANF()) GO TO 261
        GO TO 300
      ENDIF
C
C          Should never fall through
C
      GO TO 9998
C----------------------------------------------------------------------
C          Swap particles and antiparticles if IDENT(IP)<0
C          Note forced modes for antiparticles are conjugated in table.
C----------------------------------------------------------------------
300   CONTINUE
      IF(IDENT(IP).LT.0.AND.IDENT(IP).NE.-20) THEN
        DO 310 I=1,NADD
          ID1=IDENT(NPTCL+I)
          IDENT(NPTCL+I)=IDANTI(ID1)
310     CONTINUE
      ENDIF
C
C          Set IORIG and IDCAY.
C
      NPTCL=NPTCL+NADD
      IDCAY(IP)=IPACK*NSTART+NPTCL
      JETIP=IABS(IORIG(IP))/IPACK
      DO 320 I=NSTART,NPTCL
        IORIG(I)=IP
        IDCAY(I)=0
320   CONTINUE
C
C          Evolve and hadronize partons. If it fails, start over.
C
      IF(IDABS(NADD).LT.10.OR.MOD(IDENT(NPTCL),100).EQ.0) THEN
        IF(.NOT.DECJET(IP,NADD,IDABS,PREST,WDECAY,BETA,GAMMA))
     $  GO TO 1
      ENDIF
C
      RETURN
C----------------------------------------------------------------------
C          Error messages.
C----------------------------------------------------------------------
9999  CALL PRTEVT(0)
      WRITE(ITLIS,99990) NPTCL
99990 FORMAT(//5X,'ERROR IN DECAY...NPTCL > ',I6)
      RETURN
9998  CALL PRTEVT(0)
      WRITE(ITLIS,99980) IP
99980 FORMAT(//5X,'ERROR IN DECAY...NO DECAY FOUND FOR PARTICLE',I6)
      RETURN
      END
