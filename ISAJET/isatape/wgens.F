#include "isajet/pilot.h"
      SUBROUTINE WGENS
C
C          Copy event information into ZEVEL and call BUFOUT.
C          If number of words required exceeds MAXLEN-8, the number
C          of records written=no. of words/(MAXLEN-8)+1
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/mbgen.inc"
#include "isajet/keys.inc"
#include "isajet/idrun.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetset.inc"
#include "isajet/jetsig.inc"
#include "isajet/partcl.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/primar.inc"
#include "isajet/zevel.inc"
#include "isajet/final.inc"
#include "isajet/totals.inc"
#include "isajet/wsig.inc"
C
      INTEGER I2,I1,JET,K,IEX,IL,ITA,I,NPSTA
C
      ITA=IABS(ITEVT)
      IZEVEL(1)=100
      IZEVEL(2)=1
      IL=3
      CALL MOVLEI(IDVER,IZEVEL(IL),4)
      IL=IL+4
      CALL MOVLEL(KEYS(1),LZEVEL(IL),MXKEYS)
      IL=IL+MXKEYS
      IZEVEL(IL)=NJET
      IL=IL+1
      CALL MOVLEV(P(1),ZEVEL(IL),59)
      IL=IL+59
      CALL MOVLEV(SIGF,ZEVEL(IL),1)
      IL=IL+1
      IF(.NOT.KEYS(4)) THEN
        ZEVEL(IL)=SIGMA
        ZEVEL(IL+1)=SIGEVT
        ZEVEL(IL+2)=WT
        IL=IL+3
      ENDIF
C          IF ITEVT.LT.0 WRITE ONLY STABLE PARTICLES AND FLAG
C          BY NPTCL=-(NO. OF STABLE PARTICLES)
      IF(ITEVT.GT.0) THEN
        IZEVEL(IL)=NPTCL
      ELSE
        NPSTA=0
        DO 990 I=1,NPTCL
990     IF(IDCAY(I).EQ.0) NPSTA=NPSTA+1
        IZEVEL(IL)=-NPSTA
      ENDIF
      IL=IL+1
      IF(NJET.GT.0) THEN
        IEX=NJET*5
        CALL MOVLEV(PJETS(1,1),ZEVEL(IL),IEX)
        IL=IL+IEX
        CALL MOVLEI(IDJETS(1),IZEVEL(IL),NJET)
        IL=IL+NJET
      ENDIF
      IF(KEYS(3).OR.KEYS(7).OR.KEYS(11)) THEN
        CALL MOVLEV(QWJET(1),ZEVEL(IL),6)
        IL=IL+6
        CALL MOVLEV(QMW,ZEVEL(IL),16)
        IL=IL+16
        I1=JWTYP
        I2=JETTYP(3)
        IZEVEL(IL)=SIGLLQ
        IL=IL+1
      ENDIF
      IF(KEYS(6).OR.KEYS(7)) THEN
        IZEVEL(IL)=NPAIR
        IL=IL+1
        IF(NPAIR.NE.0) THEN
          CALL MOVLEV(PPAIR(1,1),ZEVEL(IL),5*NPAIR)
          IL=IL+5*NPAIR
          CALL MOVLEI(IDPAIR(1),IZEVEL(IL),NPAIR)
          IL=IL+NPAIR
          CALL MOVLEI(JPAIR(1),IZEVEL(IL),NPAIR)
          IL=IL+NPAIR
        ENDIF
      ENDIF
      IZEVEL(IL)=NJSET
      IL=IL+1
      CALL MOVLEI(NKINPT,IZEVEL(IL),5)
      IL=IL+5
      CALL MOVLEI(NPOM,IZEVEL(IL),1)
      IL=IL+1
C
C          /JETSET/ COMMON BLOCK
      IF(NJSET.LT.1) GOTO 12
      DO 50 I=1,NJSET
      CALL MOVLEV(PJSET(1,I),ZEVEL(IL),5)
      IL=IL+5
      IZEVEL(IL)=JORIG(I)
      IZEVEL(IL+1)=JTYPE(I)
      IZEVEL(IL+2)=JDCAY(I)
      IL=IL+3
      IF(IL.LE.MAXLEN-9) GO TO 50
      IZEVEL(1)=IZEVEL(1)+1
      CALL BUFOUT(IL)
      IF(I.EQ.NJSET) GO TO 12
50    CONTINUE
C
C          /PARTCL/ COMMON BLOCK
C          IF ITEVT.LT.0, WRITE OUT ONLY STABLE PARTICLES
C          FLAG BY NPTCL=-(NO. OF STABLE PARTICLES)
C          SUPPRESS ORIGIN AND DECAY INFORMATION
12    IF(NPTCL.EQ.0) GOTO 999
      IF(ITEVT.GT.0) GOTO 997
C          ONLY STABLE PARTICLES
      DO 992 K=1,NPTCL
      IF(IDCAY(K).NE.0) GOTO 992
      JET=IABS(IORIG(K))/1000
      CALL MOVLEV(PPTCL(1,K),ZEVEL(IL),5)
      IZEVEL(IL+5)=(JET*10000+IABS(IDENT(K)))*ISIGN(1,IDENT(K))
      IL=IL+6
      IF(IL.LE.MAXLEN-6) GOTO 992
      IZEVEL(1)=IZEVEL(1)+1
      CALL BUFOUT(IL)
      IF(K.EQ.NPTCL) RETURN
  992 CONTINUE
      GOTO 999
  997 CONTINUE
C          ALL PARTICLES
C          NOTE IDCAY CAN EXCEED 2**24 LIMIT OF PAIRPAK
      DO 998 K=1,NPTCL
      CALL MOVLEV(PPTCL(1,K),ZEVEL(IL),5)
      IZEVEL(IL+5)=IORIG(K)
      IZEVEL(IL+6)=IDENT(K)
      IZEVEL(IL+7)=IDCAY(K)/IPACK
      IZEVEL(IL+8)=MOD(IDCAY(K),IPACK)
      IL=IL+9
      IF(IL.LE.MAXLEN-9) GOTO 998
      IZEVEL(1)=IZEVEL(1)+1
      CALL BUFOUT(IL)
      IF(K.EQ.NPTCL) RETURN
  998 CONTINUE
  999 CONTINUE
      IZEVEL(1)=IZEVEL(1)+1
      CALL BUFOUT(IL)
      RETURN
      END
