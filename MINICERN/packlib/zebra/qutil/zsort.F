*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:44  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZSORT (IXSTOR,LGOP,JWORD)

C-    SORT BANKS AT LGO FOR WORDS Q(L+JWORD) TO BE IN INCREASING ORDER

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    JWORD(9), LGOP(9)
      REAL         KEYN, KEYPL, KEYNG
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZSOR, 4HT    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZSORT  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZSORT   ')
#endif

#include "zebra/q_jbyt.inc"

      LGO = LGOP(1)
      IF (LGO.EQ.0)          RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstore.inc"

      JW    = JWORD(1)
      KGONG = LQSTA(KQT+2) - 1
      LLNG  = KGONG
      KGOPL = KGONG - 1
      LLPL  = KGOPL
      KEYPL = -1.
      KEYNG =  0.
      IFL   =  0

      LN    = LGO
      KEYN  =  Q(KQS+LN+JW)
      IF (KEYN.LT.0.)              GO TO 62

C--------          +VE IN-SEQUENCE LOOP FOR BANKS WITH +VE KEY
C--                          KEEP GOING FOR INCREASING KEYS

   21 LQ(KQS+LLPL) = LN
      GO TO 24

   23 KEYPL = KEYN
      LLPL  = LN
      LN    = LQ(KQS+LN)
      IF (LN.EQ.0)                 GO TO 81
      KEYN  =  Q(KQS+LN+JW)
      IF (KEYN.LT.0.)              GO TO 61
   24 IF (KEYN-KEYPL)        28, 23, 23

   28 LS = KGOPL

C--------          OUT-SEQUENCE LOOP, FIND PLACE FOR BANK IN THE CHAIN
C--                          OF BANKS ALREADY SORTED, +VE OR -VE CHAIN

   41 LNX = LQ(KQS+LN)
      IFL = 7
   43 KS  = LS
      LS  = LQ(KQS+LS)
      IF (KEYN- Q(KQS+LS+JW))   48, 43, 43

   48 LQ(KQS+LN) = LS
      LQ(KQS+KS) = LN
      IF (LNX.EQ.0)                GO TO 81
      LN   = LNX
      KEYN =  Q(KQS+LN+JW)
      IF (KEYN.GE.0.)              GO TO 21

C--------          -VE IN-SEQUENCE LOOP FOR BANKS WITH -VE KEY
C--                          KEEP GOING FOR INCREASING KEYS

   61 IFL = 7
   62 LS = KGONG
      LQ(KQS+LLNG) = LN
      IF (KEYNG.NE.0.)             GO TO 64

   63 KEYNG = KEYN
      LLNG  = LN
      LN    = LQ(KQS+LN)
      IF (LN.EQ.0)                 GO TO 81
      KEYN  =  Q(KQS+LN+JW)
      IF (KEYN.GE.0.)              GO TO 21
   64 IF (KEYN-KEYNG)        41, 63, 63

C----              FINISHED, LINK +VE AND -VE STREAMS, CHAIN K-LINKS

   81 IF (IFL.EQ.0)                GO TO 999
      LQ(KQS+LLPL) = 0
      LQ(KQS+LLNG) = LQ(KQS+KGOPL)

      K = LQ(KQS+LGO+2)
      L = LQ(KQS+KGONG)
      LQ(KQS+L+2) = K
      IF (K.NE.0)  LQ(KQS+K)=L
      LGOP(1) = L

   84 K = L
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 999
      LQ(KQS+L+2) = K
      GO TO 84

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
