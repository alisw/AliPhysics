*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:26:30  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/04/18 16:11:05  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZIN (CHPA1,IPA2,IPA3,IPA4)

C-    Processor down transfer

C-     CHPA1  processor ID in A4
C-      IPA2  = 0  no further down transfer
C-            = 1 yes further down transfer
C-      IPA3  NAN = number of processor constants
C-      IPA4  extra features
C-   IPA4(2)  NCR = number of conditions to be recorded
C-   IPA4(3)  NLS = number of wsp links to be saved
C-   IPA4(4)  NDS = number of wsp data words to be saved

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      DIMENSION    IPA2(7),IPA3(7),IPA4(7)
      CHARACTER    CHPA1*4
#if defined(CERNLIB_A4)
      CHARACTER    CHIAM*4
#endif
#if defined(CERNLIB_A8)
      CHARACTER    CHIAM*8
#endif
#if defined(CERNLIB_EQUHOLCH)
      EQUIVALENCE (CHIAM, IAMID)
#endif

      DIMENSION    MMJZFO(5)
      PARAMETER   (MXREOD = 2097152)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZIN, 4H     /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZIN   /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZIN    ')
#endif
      DATA  MMJZFO / 4HJZFO, 0, 0, 1, 2 /

#include "zebra/q_jbit.inc"

#include "zebra/qtraceq.inc"
#include "zebra/qstorjz.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      CHIAM  = CHPA1
#if !defined(CERNLIB_EQUHOLCH)
      CALL UCTOH (CHIAM, IAMID,4,4)
#endif
      IAFLDW = IPA2(1)
      MINIT  = 7
      IF (LQSV.EQ.0)               GO TO 21

C----              Remember present state

      IQ(KQS+LQJZ+2*JQLEV+2) = NQLINK
      IQ(KQS+LQJZ+2*JQLEV+3) = LQSTA(KQT+1)

      N = IQ(KQS+LQSV+3)
      IF (N.NE.0) CALL UCOPY (LQ(KQS+NQREF+1),LQ(KQS+LQSV-N-3),N)
      N = IQ(KQS+LQSV+4)
      IF (N.NE.0)  THEN
          L = LQAN + IQ(KQS+LQAN)
          CALL UCOPY (LQ(KQS+NQLINK+1),IQ(KQS+L+1),N)
        ENDIF

#if defined(CERNLIB_JZTIME)
#include "zebra/jztimin.inc"
#endif

C----              Find SV bank

   21 L = LQ(KQS+LQJZ-4)
      IF (L.EQ.0)                  GO TO 24
      J = IUCOMP (IAMID,IQ(KQS+L+2),IQ(KQS+L+1))
      IF (J.EQ.0)                  GO TO 24
      LQSV = LQ(KQS+L-J)
      GO TO 25

   24 LQSV   = LZFIND (IXSTJZ,LQ(KQS+LQJZ-3), IAMID,1)
      IF (LQSV.EQ.0)               GO TO 81
   25 IQ(KQS+LQSV+2) = IQ(KQS+LQSV+2) + 1
      LCD  = LQSV + JQNACC
      LQAN = LCD + IQ(KQS+LCD) + 1

C--                Copy flags

#if defined(CERNLIB_QDEBUG)
      CALL VZERO (JQFLAG,JQMFLW)
      L = LQAN + IQ(KQS+LQAN) + 1 + IQ(KQS+LQSV+4)
      N = IQ(KQS+L)
      IF (N.NE.0)  CALL UCOPY (IQ(KQS+L+1),JQFLAG,N)

      IF (JBIT(JQLLEV,9).NE.0)  WRITE (IQLOG,9024) JQLEV,NQME(1),
     +                          IQ(KQS+LQSV+1),IQ(KQS+LQSV+2)

 9024 FORMAT (/' =======  JZIN   level',I2,', "',A4,' down to "',A4,
     FI8,'th entry')
#endif

      JQEALL = JQEALL + 1
      IF (JQEALL.EQ.0)             GO TO 71

C----              Step level

   31 NQME(1) = IQ(KQS+LQSV+1)
      JQLEV   = JQLEV + 1
      IF (JQLEV.GT.JQMLEV)         GO TO 91
      J     = LQJZ - JQLEV - 6
      LQDW  = LQ(KQS+J)
      LQUP  = LQ(KQS+J+1)
      J     = J - JQMLEV
      LQ(KQS+J) = LQSV
      IQUEST(1) = MINIT
      IF (IAFLDW.NE.0)             GO TO 37
      LQDW = 0
#include "zebra/qtrace99.inc"
      RETURN

C--                Clear the down call bank

   37 IF (JQLEV.EQ.JQMLEV)         GO TO 92
      CALL VZERO (LQ(KQS+LQDW-JQCBNL),JQCBNL)
      CALL VZERO (IQ(KQS+LQDW+1),     JQCBND)
      IQ(KQS+LQDW) = MSBYT (0, IQ(KQS+LQDW),1,18)
      GO TO 999

C----              Re-order SV structure every now and then

   71 JQREOD = MIN (4*JQREOD,MXREOD)
      JQEALL = -JQREOD
      IF (JQREOD.GE.MXREOD)        GO TO 31
      L = LQ(KQS+LQJZ-3)
      CALL ZTOPSY (IXSTJZ,L)
      CALL ZSORTI (IXSTJZ,L,2)
      CALL ZTOPSY (IXSTJZ,L)
      NPR  = NZBANK (IXSTJZ,L)
      INC  = 10

      LFO = LQ(KQS+LQJZ-4)
      IF (LFO.EQ.0)                GO TO 72
      IF (NPR.LE.IQ(KQS+LFO-1))       GO TO 74
      CALL MZDROP (IXSTJZ,IQ(KQS+LFO), '.')
      INC = 4

   72 MMJZFO(2) = NPR + INC
      MMJZFO(4) = MMJZFO(2) + 1
      CALL MZLIFT (IXDVJZ,LFO,LQJZ,-4,MMJZFO,0)
      JQREOD = 512
      JQEALL = -JQREOD

   74 IQ(KQS+LFO+1) = NPR
      J = 0
      L = LQJZ - 3
   75 L = LQ(KQS+L)
      IF (L.EQ.0)                  GO TO 31
      J = J + 1
      LQ(KQS+LFO-J)   = L
      IQ(KQS+LFO+J+1) = IQ(KQS+L+1)
      GO TO 75

C----              Processor not yet initialized

   81 IANAN = IPA3(1)
      IANCR = 10
      IANLSV = 0
      IANDSV = 0

      N = IPA4(1)
      IF (N.GE.0)  THEN
          N =  MIN (N,3)
          CALL UCOPY (IPA4(2),IANCR,N)
        ENDIF
      CALL JZLIFT
      MINIT = IQUEST(1)
      LQSV  = IQUEST(2)

      LFO = LQ(KQS+LQJZ-4)
      IF (LFO.EQ.0)                GO TO 25
      NFO = IQ(KQS+LFO+1) + 1
      IF (NFO.GE.IQ(KQS+LFO-1))       GO TO 25
      LQ(KQS+LFO-NFO)   = LQSV
      IQ(KQS+LFO+NFO+1) = IAMID
      IQ(KQS+LFO+1)     = NFO
      GO TO 25

C----              JQMLEV exeeded

   92 NQCASE = 1
   91 NQCASE = NQCASE + 1
      NQFATA = 3
      IQUEST(11) = IAMID
      IQUEST(12) = IAFLDW
      IQUEST(13) = JQLEV
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
