*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:27  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1996/04/18 16:13:45  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if (defined(CERNLIB_VAX))||((defined(CERNLIB_LINUX) && !defined(CERNLIB_PPC)))||(defined(CERNLIB_QMVMI))||(defined(CERNLIB_QMDOS))
#include "qutyvax/zsorth.F"
#else
      SUBROUTINE ZSORTH (IXSTOR,LGOP,JWORD)

C-    SORT BANKS AT LGO FOR WORDS Q(L+JWORD) TO BE IN INCREASING ORDER

C-    ALPHABETIC SORT DONE WITH INTEGER COMPARISON
C-    THE SIGN-BIT IS CONSIDERED AS AN ORDINARY BIT, FOR CONSISTENCY
C-    ON ANY GIVEN MACHINE THIS MAY NOT BE VERY SATISFACTORY

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    JWORD(9), LGOP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZSOR, 4HTH   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZSORTH /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZSORTH  ')
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
      KEYPL = -1
      KEYNG =  0
      IFL   =  0

      LN    = LGO
      KEYN  = IQ(KQS+LN+JW)
      IF (KEYN.GE.0 )              GO TO 22
      GO TO 62

C--------          +VE IN-SEQUENCE LOOP FOR BANKS WITH +VE KEY
C--                          KEEP GOING FOR INCREASING KEYS

   21 IFL = 7
   22 LQ(KQS+LLPL) = LN
      GO TO 24

   23 KEYPL = KEYN
      LLPL  = LN
      LN    = LQ(KQS+LN)
      IF (LN.EQ.0)                 GO TO 81
      KEYN  = IQ(KQS+LN+JW)
      IF (KEYN.LT.0 )              GO TO 62
   24 IF (KEYN-KEYPL)        28, 23, 23

   28 LS = KGOPL

C--------          OUT-SEQUENCE LOOP, FIND PLACE FOR BANK IN THE CHAIN
C--                          OF BANKS ALREADY SORTED, +VE OR -VE CHAIN

   41 LNX = LQ(KQS+LN)
      IFL = 7
   43 KS  = LS
      LS  = LQ(KQS+LS)
      IF (KEYN-IQ(KQS+LS+JW))   48, 43, 43

   48 LQ(KQS+LN) = LS
      LQ(KQS+KS) = LN
      IF (LNX.EQ.0)                GO TO 81
      LN   = LNX
      KEYN = IQ(KQS+LN+JW)
      IF (KEYN.GE.0 )              GO TO 22

C--------          -VE IN-SEQUENCE LOOP FOR BANKS WITH -VE KEY
C--                          KEEP GOING FOR INCREASING KEYS

   62 LS = KGONG
      LQ(KQS+LLNG) = LN
      IF (KEYNG.NE.0 )             GO TO 64

   63 KEYNG = KEYN
      LLNG  = LN
      LN    = LQ(KQS+LN)
      IF (LN.EQ.0)                 GO TO 81
      KEYN  = IQ(KQS+LN+JW)
      IF (KEYN.GE.0 )              GO TO 21
   64 IF (KEYN-KEYNG)        41, 63, 63

C----              FINISHED, LINK +VE AND -VE STREAMS, CHAIN K-LINKS

   81 IF (IFL.EQ.0)                GO TO 999
      LQ(KQS+LLNG) = 0
      LQ(KQS+LLPL) = LQ(KQS+KGONG)

      K = LQ(KQS+LGO+2)
      L = LQ(KQS+KGOPL)
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
#endif
