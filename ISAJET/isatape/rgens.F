#include "isajet/pilot.h"
      SUBROUTINE RGENS(IFLAG)
C
C          Inverse of WRGEN
C          Read a record by a call BUFIN
C          If record type is not event type return.
C          If RGENS called with IFLAG=10 return without unpacking.
C          Unpack ZEVEL into appropriate common blocks.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/mbgen.inc"
#include "isajet/keys.inc"
#include "isajet/rectp.inc"
#include "isajet/idrun.inc"
#include "isajet/jetpar.inc"
#include "isajet/jetset.inc"
#include "isajet/jetsig.inc"
#include "isajet/partcl.inc"
#include "isajet/pjets.inc"
#include "isajet/pinits.inc"
#include "isajet/primar.inc"
#include "isajet/zevel.inc"
#include "isajet/totals.inc"
#include "isajet/wsig.inc"
#include "isajet/final.inc"
C
      INTEGER IFLAG
      INTEGER I,IEX,IZ5,K,IFL,ISAV,IL
C
      IFL=IFLAG
      CALL BUFIN(IL,IFLAG)
      IF(IFLAG.NE.0) RETURN
      ISAV=IZEVEL(1)/100
      IRECTP=ISAV*100
      IREC=MOD(IZEVEL(1),100)
      IF(IRECTP.EQ.200) RETURN
      IF(IRECTP.EQ.300) RETURN
      IF(IFL.EQ.10) RETURN
      IL=3
      CALL MOVLEI(IZEVEL(IL),IDVER,4)
      IL=IL+4
      CALL MOVLEL(LZEVEL(IL),KEYS(1),MXKEYS)
      IL=IL+MXKEYS
      NJET=IZEVEL(IL)
      IL=IL+1
      CALL MOVLEV(ZEVEL(IL),P(1),59)
      IL=IL+59
      CALL MOVLEV(ZEVEL(IL),SIGF,1)
      IL=IL+1
      IF(.NOT.KEYS(4)) THEN
        SIGMA=ZEVEL(IL)
        SIGEVT=ZEVEL(IL+1)
        WT=ZEVEL(IL+2)
        IL=IL+3
      ENDIF
      NPTCL=IZEVEL(IL)
      IL=IL+1
      IF(NJET.GT.0) THEN
        IEX=NJET*5
        CALL MOVLEV(ZEVEL(IL),PJETS(1,1),IEX)
        IL=IL+IEX
        CALL MOVLEI(IZEVEL(IL),IDJETS(1),NJET)
        IL=IL+NJET
      ENDIF
      IF(KEYS(3).OR.KEYS(7).OR.KEYS(11)) THEN
        CALL MOVLEV(ZEVEL(IL),QWJET(1),6)
        IL=IL+6
        CALL MOVLEV(ZEVEL(IL),QMW,16)
        IL=IL+16
        SIGLLQ=ZEVEL(IL)
        IL=IL+1
      ENDIF
      IF(KEYS(6).OR.KEYS(7)) THEN
        NPAIR=IZEVEL(IL)
        IL=IL+1
        IF(NPAIR.NE.0) THEN
          CALL MOVLEV(ZEVEL(IL),PPAIR(1,1),5*NPAIR)
          IL=IL+5*NPAIR
          CALL MOVLEI(IZEVEL(IL),IDPAIR(1),NPAIR)
          IL=IL+NPAIR
          CALL MOVLEI(IZEVEL(IL),JPAIR(1),NPAIR)
          IL=IL+NPAIR
        ENDIF
      ENDIF
      NJSET=IZEVEL(IL)
      IL=IL+1
      CALL MOVLEI(IZEVEL(IL),NKINPT,5)
      IL=IL+5
      CALL MOVLEI(IZEVEL(IL),NPOM,1)
      IL=IL+1
C
C          /JETSET/ COMMON BLOCK
      IF(NJSET.LT.1) GOTO 12
      DO 50 I=1,NJSET
      CALL MOVLEV(ZEVEL(IL),PJSET(1,I),5)
      IL=IL+5
      JORIG(I)=IZEVEL(IL)
      JTYPE(I)=IZEVEL(IL+1)
      JDCAY(I)=IZEVEL(IL+2)
      IL=IL+3
      IF(IL.LE.MAXLEN-9) GO TO 50
      IF(I.EQ.NJSET.AND.NPTCL.EQ.0) GO TO 12
      CALL BUFIN(IL,IFLAG)
      IF(IFLAG.NE.0) RETURN
50    CONTINUE
C
C          /PARTCL/ COMMON BLOCK
C          NPTCL.LT.0 IMPLIES ONLY STABLE PARTICLES ON THIS FILE
C          ORIGIN AND DECAY INFORMATION SUPPRESSED
12    IF(NPTCL.EQ.0) GOTO 999
      IF(NPTCL.GT.0) GOTO 997
C          ONLY STABLE PARTICLES
      NPTCL=-NPTCL
      DO 992 K=1,NPTCL
      CALL MOVLEV(ZEVEL(IL),PPTCL(1,K),5)
      IZ5=IABS(IZEVEL(IL+5))
      IORIG(K)=(IZ5/10000)*1000
      IDENT(K)=MOD(IZ5,10000)*ISIGN(1,IZEVEL(IL+5))
      IDCAY(K)=0
      IL=IL+6
      IF(IL.LE.MAXLEN-6) GOTO 992
      IF(K.EQ.NPTCL) RETURN
      CALL BUFIN(IL,IFLAG)
      IF(IFLAG.NE.0) RETURN
  992 CONTINUE
      RETURN
C          ALL PARTICLES
C          NOTE THAT IDCAY CAN EXCEED 2**24 LIMIT OF PAIRPAK
  997 CONTINUE
      DO 998 K=1,NPTCL
      CALL MOVLEV(ZEVEL(IL),PPTCL(1,K),5)
      IORIG(K)=IZEVEL(IL+5)
      IDENT(K)=IZEVEL(IL+6)
      IDCAY(K)=IZEVEL(IL+7)*IPACK+IZEVEL(IL+8)
      IL=IL+9
      IF(IL.LE.MAXLEN-9) GOTO 998
      IF(K.EQ.NPTCL) RETURN
      CALL BUFIN(IL,IFLAG)
      IF(IFLAG.NE.0) RETURN
  998 CONTINUE
  999 CONTINUE
      RETURN
      END
