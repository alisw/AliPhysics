*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:26:32  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/04/18 16:12:37  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZRELB

C-    Relocator for links in banks

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZRE, 4HLB   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZRELB /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZRELB  ')
#endif

#if !defined(CERNLIB_QDEBUG)
#include "zebra/q_jbyt.inc"
#endif


#include "zebra/qtrace.inc"

      LFIXLO = LQ(LQTA-1)
      LFIXRE = LQ(LQTA)
      LFIXHI = LQ(LQTE)
      JHIGO  = (LQTE-LQTA) / 4
      NENTR  = JHIGO - 1

      IF (NENTR.EQ.0)  THEN
          LADTB1 = LQ(LQTA+1)
          NRLTB2 = LQ(LQTA+2)
          IFLTB3 = LQ(LQTA+3)
        ENDIF

      LMRNX  = LQMTA
   12 LMR    = LMRNX
      IF (LMR.GE.LQMTE)            GO TO  999
      LMRNX  = LMRNX + 8

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9813) (LQ(J+LMR-1),J=1,8)
 9813 FORMAT (1X/' DEVZE MZRELB.  DIV, ACT,   NSH, LF, LL, LTA, LTB'
     F,', NFREE'/14X,2I5,6I7)
#endif

      IACT   = LQ(LMR+1)
      IF (IACT.LE.0)               GO TO 12
      IF (IACT.EQ.4)               GO TO 12
      LSTOP  = LQ(LMR+4)
      IF (IACT.EQ.3)               GO TO 14
      LN     = LQ(LMR+3)
      LDEAD  = LSTOP
      GO TO 19

   14 LSEC   = LQRTA + LQ(LMR+5) - 4

C------            Next bank,  check if dead group

   16 LSEC   = LSEC + 4
      LNX    = LQ(LSEC)
      LDEAD  = LQ(LSEC+1)

   17 LN = LNX
      IF (LN.GE.LSTOP)             GO TO 12
      IF (LN.EQ.LDEAD)             GO TO 16

C--                Next bank,  alive

   19 CONTINUE

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLN (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 91
      LNX = IQNX
      IF (IQND.LT.0)               GO TO 17
      LS  = IQLS
      LX  = LS + 3
      L2  = LS - IQNS
      L1  = LS - IQNL
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9819) LN, IQID,IQNL,IQNS,IQND
 9819 FORMAT (1X/' DEVZE MZRELB.  Do bank at LN =',I7
     F,'  ID,NL,NS,ND= ',A4,3I7)
#endif
#if !defined(CERNLIB_QDEBUG)
      NST = JBYT (LQ(KQS+LN),1,16) - 11
      IF (NST.LT.0)  THEN
          LNX = LN + NST + 11
          GO TO 17
        ELSE
          LS  = LN + NST
          LX  = LS + 3
          L2  = LS - IQ(KQS+LS-2)
          L1  = LS - IQ(KQS+LS-3)
          LNX = LS + IQ(KQS+LS-1) + 9
        ENDIF
#endif
      IF   (NENTR)           66, 46, 26

C--------------    2 OR MORE RELOCATION INTERVALS       -------------

C----              Next link

   24 LQ(KQS+L1)= 0

   25 L1 = L1 + 1
      IF (L1.EQ.LX)                GO TO 17
   26 LFIRST= LQ(KQS+L1)
   27 LINK  = LQ(KQS+L1)
      IF (LINK.EQ.0)               GO TO 25

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9827) LINK,L1
 9827 FORMAT (16X,'Link =',I7,' from L1 =',I7)
#endif
      IF (IQFLIO.EQ.0)  THEN
          IF (LINK.LT.LFIXLO)      GO TO 25
          IF (LINK.GE.LFIXHI)      GO TO 25
          IF (LINK.LT.LFIXRE)      GO TO 24
        ELSE
          IF (LINK.LT.LFIXRE)      GO TO 24
          IF (LINK.GE.LFIXHI)      GO TO 24
        ENDIF

C--                Binary search in relocator table

      JLOW = 0
      JHI  = JHIGO

   29 JEX = (JHI+JLOW) / 2
      IF (JEX.EQ.JLOW)             GO TO 31
      IF (LINK.GE.LQ(LQTA+4*JEX))  GO TO 30
      JHI  = JEX
      GO TO 29

   30 JLOW = JEX
      GO TO 29

C--                Relocate

   31 JTB = LQTA + 4*JLOW

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9831) JLOW, (LQ(JTB+J-1),J=1,4)
 9831 FORMAT (50X,'Entry',I5,',',4I7)
#endif
      IF (LINK.GE.LQ(JTB+1))             GO TO 33
      LQ(KQS+L1) = LINK + LQ(JTB+2)
      GO TO 25

C----              Link into dead area

   33 IF   (LQ(JTB+3))       25, 24, 34

C--                Bridge structural link

   34 IF (L1.LT.L2)                GO TO 24
      IF (LS+1-L1)           36, 24, 35
   35 CONTINUE

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.NE.0)             GO TO 92
#endif
      LINK = LQ(KQS+LINK)
      LQ(KQS+L1) = LINK
      IF (LINK.NE.LFIRST)          GO TO 27
      GO TO 24

C--                Reverse bridging of s-link

   36 LINK = LQ(KQS+LINK+2)
      LQ(KQS+L1) = LINK
      GO TO 27

C--------------    1 RELOCATION INTERVAL ONLY           -------------

C----              Next link

   44 LQ(KQS+L1)= 0

   45 L1 = L1 + 1
      IF (L1.EQ.LX)                GO TO 17
   46 LFIRST= LQ(KQS+L1)
   47 LINK  = LQ(KQS+L1)
      IF (LINK.EQ.0)               GO TO 45

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9827) LINK,L1
#endif
      IF (IQFLIO.EQ.0)  THEN
          IF (LINK.LT.LFIXLO)      GO TO 45
          IF (LINK.GE.LFIXHI)      GO TO 45
          IF (LINK.LT.LFIXRE)      GO TO 44
          IF (LINK.GE.LADTB1)      GO TO 53
        ELSE
          IF (LINK.LT.LFIXRE)      GO TO 44
          IF (LINK.GE.LADTB1)      GO TO 44
        ENDIF

C--                Relocate

      LQ(KQS+L1) = LINK + NRLTB2
      GO TO 45

C----              Link into dead area

   53 IF   (IFLTB3)           45, 44, 54

C--                Bridge structural link

   54 IF (L1.LT.L2)                GO TO 44
      IF (LS+1-L1)           56, 44, 55
   55 CONTINUE

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.NE.0)             GO TO 92
#endif
      LINK = LQ(KQS+LINK)
      LQ(KQS+L1) = LINK
      IF (LINK.NE.LFIRST)          GO TO 47
      GO TO 44

C--                Reverse bridging of s-link

   56 LINK = LQ(KQS+LINK+2)
      LQ(KQS+L1) = LINK
      GO TO 47

C--------------    NO RELOCATION INTERVAL               -------------

C----              Next link

   64 LQ(KQS+L1)= 0

   65 L1 = L1 + 1
      IF (L1.EQ.LX)                GO TO 17
   66 LINK  = LQ(KQS+L1)
      IF (LINK.EQ.0)               GO TO 65

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)   WRITE (IQLOG,9827) LINK,L1
#endif
      IF (LINK.LT.LFIXLO)          GO TO 65
      IF (LINK.GE.LFIXHI)          GO TO 65
      GO TO 64

C------            Error conditions
#if defined(CERNLIB_QDEBUG)
   92 NQCASE = 1
      NQFATA = 2
      LN     = LS
      IQUEST(12) = L1
      IQUEST(13) = LINK
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(11) = LN
      IF (IQFLIO.NE.0)             GO TO 98
#include "zebra/qtofatal.inc"

   98 IQUEST(9) = NQCASE
      IQUEST(10)= NQFATA
      NQCASE = 0
      NQFATA = 0
      IQFLIO = -7
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
