#include "isajet/pilot.h"
      SUBROUTINE ISAWBG
C
C          Write initial record (type 200)
C          Inverse of RDBEG
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/dylim.inc"
#include "isajet/frgpar.inc"
#include "isajet/idrun.inc"
#include "isajet/jetlim.inc"
#include "isajet/keys.inc"
#include "isajet/primar.inc"
#include "isajet/qcdpar.inc"
#include "isajet/qlmass.inc"
#include "isajet/q1q2.inc"
#include "isajet/types.inc"
#include "isajet/xmssm.inc"
C
#include "isajet/zvout.inc"
#include "isajet/zevel.inc"
#include "isajet/final.inc"
C
      INTEGER NL,IL,ITA
C
C       Keep entry point WRBEG for backward compatibility
      ENTRY WRBEG
C
      ITA=IABS(ITEVT)
      NKINF=0
      SIGF=0.
      ALUM=0.
      ACCEPT=0.
      NRECS=0
#if defined(CERNLIB_CDCPACK)
      CALL ZEROL(ZVOUT,512)
#endif
      CALL ZEROL(ZEVEL,MAXLEN)
      IL=3
      CALL MOVLEI(IDVER,IZEVEL(IL),4)
      IL=IL+4
      CALL MOVLEI(NJET,IZEVEL(IL),7)
      IL=IL+7
      NL=NJET*MXGOQ
      IF(NJET.NE.0) CALL MOVLEL(GOQ(1,1),LZEVEL(IL),NL)
      IL=IL+NL
      CALL MOVLEL(KEYS(1),LZEVEL(IL),10)
      IL=IL+10
      CALL MOVLEV(PMIN(1),ZEVEL(IL),36)
      IL=IL+36
      IF(.NOT.KEYS(3)) GOTO 11
      CALL MOVLEV(QMIN,ZEVEL(IL),12)
      IL=IL+12
   11 CONTINUE
      CALL MOVLEL(GODY(1),LZEVEL(IL),5)
      IL=IL+5
      CALL MOVLEV(PUD,ZEVEL(IL),22)
      IL=IL+22
      CALL MOVLEV(ALAM,ZEVEL(IL),4)
      IL=IL+4
      CALL MOVLEV(AMLEP(6),ZEVEL(IL),3)
      IL=IL+3
      CALL MOVLEI(LOC(1),IZEVEL(IL),100)
      IL=IL+100
      CALL MOVLEL(GOMSSM,LZEVEL(IL),1)
      IL=IL+1
      CALL MOVLEV(XGLSS,ZEVEL(IL),11)
      IL=IL+11
      CALL MOVLEL(GOSUG,LZEVEL(IL),1)
      IL=IL+1
      CALL MOVLEV(XM0SU,ZEVEL(IL),5)
      IL=IL+5
C
      IZEVEL(1)=200
      IZEVEL(2)=1
      CALL BUFOUT(IL)
      RETURN
      END
