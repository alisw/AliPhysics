*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:26:27  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/04/18 16:10:48  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOFFX (IUHEAD)

C-    Write to buffer operations for file format Exchange,
C-    subsidiary to FZOUT

C-    Controlling parameter : IDX(2)
C-
C-    IDX(2)   = 1   write start/end-of-run
C-             > 1   write pilot for d/s
C-             = 0   write bank material for d/s
C-             = -1  flush the buffer
C-             = -2  End-of-File
C-             = -3  End-of-Data

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzioc.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcseg.inc"
C--------------    End CDE                             --------------
      DIMENSION    IUHEAD(99)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOF, 4HFX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOFFX  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOFFX  ')
#endif
#if !defined(CERNLIB_FQXISN)
#include "zebra/q_jbyt.inc"
#endif

#include "zebra/qtraceq.inc"

      IF   (IDX(2))          71, 41, 11
   11 IF (IDX(2).EQ.1)             GO TO 61

C-----------------------------------------------------------
C------            WRITE PILOT INFORMATION, STARTING LR AND D/S
C-----------------------------------------------------------

      IDX(1) = 10 + NWUHCX + NWSEGX + NWTXX + NWTABX + NWBKX
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOX.EQ.0)             GO TO 31

C-----------------------------------------------------------
C--                Pilot in exchange data format
C-----------------------------------------------------------

      MFO(1) = 3
      MFO(2) = 1
      MFO(3) = 1
      MFO(4) = -1
      JFOEND = 4
      CALL FZOTRX (IPILX, 10+NWIOX)

C--                User header

      IF (NWUHCX.NE.0)  THEN
          CALL MZIOCR (IOCHX)
          CALL FZOTRX (IUHEAD, NWUHX)
        ENDIF

C--                Segment table

      IF (NWSEGX.NE.0)  THEN
          MFO(1) = 5
          MFO(2) = -1
          JFOEND = 2
          CALL FZOTRX (IQSEGH, 2*NQSEG)
          MFO(1) = 1
          MFO(2) = -1
          CALL FZOTRX (IQSEGD, NQSEG)
        ENDIF

C--                Text vector

      IF (NWTXX.NE.0)  THEN
          MFO(1) = 5
          MFO(2) = 0
          JFOEND = 2
          JFOREP = 0
          CALL FZOTRX (IQ(KQSP+LTEXTX+5), NWTXX)
        ENDIF

C--                Relocation table, only FZOUT, not FZCOPY

      IF (ICOPYX.NE.0)             GO TO 39
      IF (NWTABX.NE.0)  THEN
          MFO(1) = 1
          MFO(2) = -1
          JFOEND = 2
          CALL FZOTRX (LQ(LQTA), NWTABX)
        ENDIF
      GO TO 39
#endif

C-----------------------------------------------------------
C--                Pilot in native data format
C-----------------------------------------------------------

   31 CALL FZOTRN (IPILX,10+NWIOX)

C--                User header

      IF (NWUHCX.NE.0)  CALL FZOTRN (IUHEAD, NWUHX)

C--                Segment table

      IF (NWSEGX.NE.0)  THEN
          CALL FZOTRN (IQSEGH, 2*NQSEG)
          CALL FZOTRN (IQSEGD,   NQSEG)
        ENDIF

C--                Text vector

      IF (NWTXX.NE.0)  THEN
          CALL FZOTRN (IQ(KQSP+LTEXTX+5), NWTXX)
        ENDIF

C--                Relocation table, only FZOUT, not FZCOPY

      IF (ICOPYX.NE.0)             GO TO 39
      IF (NWTABX.NE.0)  CALL FZOTRN (LQ(LQTA), NWTABX)
   39 IF (NWBKX.NE.0)              GO TO 999
      GO TO 991

C-----------------------------------------------------------
C--                WRITE BANK MATERIAL
C-----------------------------------------------------------

   41 LTB    = LQTA
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOX.EQ.0)             GO TO 51

C-----------------------------------------------------------
C--                Bank groups in exchange data format
C-----------------------------------------------------------

   42 L   = LQ(LTB)
      LE  = LQ(LTB+1)

C--                Next bank

   44 IQLN = L
      IWD  = LQ(KQS+L)
      NST  = JBYT (IWD,1,16) - 12
      IF (NST.LT.0)                GO TO 48

C--                True bank

      IQLS  = L + NST + 1
      IQNIO = JBYT (IQ(KQS+IQLS),19,4)
      IQNL  = IQ(KQS+IQLS-3)
      IQND  = IQ(KQS+IQLS-1)
      IQNX  = IQLS + IQND + 9

C-       first word, I/O words, links, next-link, up-link
      N = IQNIO + IQNL
      MFO(1) = 1
      MFO(2) = N + 3

C-       origin, numeric ID
      MFO(3) = 2
      MFO(4) = 2

C-       Hollerith ID
      MFO(5) = 5
      MFO(6) = 1

C-       NL, NS, ND, status
      MFO(7) = 1
      MFO(8) = -1
      JFOEND = 8
      CALL FZOTRX (LQ(KQS+L), N+10)

C-       data words
      IF (IQND.EQ.0)               GO TO 46
      CALL MZIOCR (LQ(KQS+IQLN))
      CALL FZOTRX (IQ(KQS+IQLS+1), IQND)

   46 L = IQNX
      IF (L.LT.LE)                 GO TO 44
      LTB = LTB + 2
      IF (LTB.LT.LQTE)             GO TO 42
      GO TO 991

C--                Short dead region

   48 NWD  = NST + 12
      IQNX = L + NWD
      MFO(1) = 1
      MFO(2) = -1
      JFOEND = 2
      CALL FZOTRX (LQ(KQS+L), NWD)
      GO TO 46
#endif

C-----------------------------------------------------------
C--                Bank-groups in native data format
C-----------------------------------------------------------

   51 L   = LQ(LTB)
      LE  = LQ(LTB+1)
      CALL FZOTRN (LQ(KQS+L), LE-L)
      LTB = LTB + 2
      IF (LTB.LT.LQTE)             GO TO 51

  991 CALL FZOREC
#include "zebra/qtrace99.inc"
      RETURN

C-----------------------------------------------------------
C--                WRITE START/END-OF-RUN
C-----------------------------------------------------------

   61 JRUN  = IQUEST(11)
      NWUHU = IDX(1) - 1
#if !defined(CERNLIB_FQXISN)
      IF (IDAFOX.EQ.0)             GO TO 64
      MFO(1) = 2
      MFO(2) = -1
      JFOEND = 2
      CALL FZOTRX (JRUN,1)
      IF (NWUHU.EQ.0)              GO TO 991
      CALL FZOTRX (IUHEAD,NWUHU)
      GO TO 991
#endif
   64 CALL FZOTRN (JRUN,1)
      IF (NWUHU.NE.0)  CALL FZOTRN (IUHEAD,NWUHU)
      GO TO 991

C-----------------------------------------------------------
C--                FLUSH BUFFER / ENDFILE
C-----------------------------------------------------------

   71 GO TO 991
      END
*      ==================================================
#include "zebra/qcardl.inc"
