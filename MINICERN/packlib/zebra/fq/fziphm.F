*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:40  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZMEMORY)
      SUBROUTINE FZIPHM

C-    Read next physical record from memory           .
C-    Service routine to FZIN, called only via FZIREC

C-    Input :      IFLAGI = 0  normal read
C-                         -1  recover to next steering block

C-                 N4SKII      is used for rapid skip of fast blocks

C-    Output :     IFLAGI = 0  all is well
C-                             otherwise : ready for re-start

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
#include "fzhci.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LRTYP,IDI(2))
*      Declaratives, DIMENSION etc.
#include "fziphrd1.inc"
* Ignoring t=pass

#include "fzstamp.inc"
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIP, 4HHM   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIPHM /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIPHM  ')
#endif

*      Declaratives, DATA
#include "fziphrd2.inc"


#include "zebra/qtrace.inc"

C----              Prepare to transfer next record

      LBUF   = IQ(KQSP+LQFI+1)
      NWMREC = IQ(KQSP+LBPARI+1)
      NFASTI = IQ(KQSP+LBPARI-5)
      NRSKIP = 0
      NBLK   = IQ(KQSP+LQFI+22)

C--                Skip complete records

      N4SKIP = MIN (N4SKII,N4RESI)
      IF (N4SKIP.LT.MAXREI)        GO TO 20
      IF (NFASTI.EQ.0)             GO TO 20
      NRSKIP = MIN (NFASTI, N4SKIP / MAXREI)
      NW4S   = NRSKIP * MAXREI
      NWMS   = NRSKIP * NWMREC
      LBUF   = LBUF   + NWMS
      N4SKII = N4SKII - NW4S
      N4RESI = N4RESI - NW4S
      NWRDAI = NWRDAI + NW4S
      NFASTI = NFASTI - NRSKIP
      NBLK   = NBLK   + NRSKIP

   20 NW4IN  = MAXREI
      NWMIN  = NWMREC

      IQ(KQSP+LBPARI-6) = NFASTI

C----              READ ONE PHYSICAL RECORD

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +      WRITE (IQLOG,9031) NBLK+1,NW4IN,NWMIN,NFASTI,NRSKIP
 9031 FORMAT (1X/' FZIPHM-  Reading Block',I7,
     F', NW32,NWmach,NRfast,NRskip=',4I6)
#endif

C--                Copy, with unpacking or byte-swop if nec.

#if defined(CERNLIB_FQNEEDCV)
#include "fziphm42.inc"
#endif
   44 CALL UCOPY (LQ(LBUF),LQ(L4STAI),NW4IN)
   46 CONTINUE

C--                Short/full dump of record read

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  CALL FZIDUM (LQ(L4STAI),NW4IN)
#endif

C-----------------------------------------------------------
C----              FAST RECORD EXPECTED
C-----------------------------------------------------------

   51 IF (NFASTI.EQ.0)             GO TO 61
      NFASTI = NFASTI - 1
      N4DONI = 0
      N4ENDI = MIN (N4RESI,MAXREI)

   58 NWRDAI = NWRDAI + MAXREI
      IQ(KQSP+LBPARI-5) = NFASTI
      IQ(KQSP+LQFI+22)  = NBLK + 1
      IQ(KQSP+LQFI+1)   = LBUF + NWMIN
      IFLAGI = 0
#include "zebra/qtrace99.inc"
      RETURN

C-----------------------------------------------------------
C----              STEERING RECORD EXPECTED
C-----------------------------------------------------------

   61 CALL FZICHH (0, LQ(L4STAI),0)
      IF (IQUEST(1).NE.0)          GO TO 64

      N4DONI = 8
      N4ENDI = NTLRI
      IF (N4ENDI.EQ.0)  N4ENDI=MAXREI
      IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      NFASTI = NFSTI
      GO TO 58

C----              Physical record length mis-match
C--                Reset if possible

   64 IF (IQUEST(1).NE.3)          GO TO 801
      IF (LBUF.NE.IQ(KQSP+LQFI+8))    GO TO 802
      CALL FZIPRL (0)
      IF (IQUEST(1).EQ.2)          GO TO 803
      IF (IQUEST(1).NE.0)          GO TO 802
      NWMREC = IQ(KQSP+LBPARI+1)
      GO TO 20


C-----------------------------------------------------------
C-                 ERROR CONDITIONS
C-----------------------------------------------------------

C-    JERROR = 521  Block header faulty
  801 JERROR = 521
      GO TO 817

C-    JERROR = 522  Block size does not match expectation
  802 JERROR = -1

C-    JERROR = 523  Block size larger than buffer
  803 JERROR = JERROR + 523
      IQUEST(14) = MAXREI
      IQUEST(15) = NWRI
      NWERR = 2

  817 IQ(KQSP+LBPARI-9) = -1
      IQ(KQSP+LBPARI-6) = 0
      IQ(KQSP+LBPARI-5) = 0
      IQ(KQSP+LBPARI-1) = 0

      JRETCD = 6
      IFLAGI = 1
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
