*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:41  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIPHR

C-    Read next physical record in sequential mode
C-    Service routine to FZIN, called only via FZIREC

C-    Input :      IFLAGI = 0  normal read
C-                         -1  recover to next steering block
C-                         -2  start

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

*      Declaratives, DIMENSION etc.
#include "fziphrd1.inc"
* Ignoring t=pass

#include "fzstamp.inc"
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIP, 4HHR   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIPHR /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIPHR  ')
#endif

*      Declaratives, DATA
#include "fziphrd2.inc"
* Ignoring t=pass


#include "zebra/qtrace.inc"

C----              DECIDE START ADR AND LENGTH FOR READ

      NWMREC = IQ(KQSP+LBPARI+1)
      NFASTI = IQ(KQSP+LBPARI-5)
      N4SKIP = MIN (N4SKII,N4RESI)
      NRSKIP = 0
      IF (N4SKIP.GE.MAXREI)  THEN
          IF (NFASTI.NE.0)  NRSKIP = MIN (NFASTI, N4SKIP / MAXREI)
        ENDIF

   20 NW4IN  = MAXREI
      NWMIN  = NWMREC

      IQ(KQSP+LBPARI-6) = NFASTI

      LIN = L4STAI
#if !defined(CERNLIB_QREADFULL)
      IF (NRSKIP.GE.2)  THEN
          NW4IN  = 90
          NWMIN  = NW4IN
        ENDIF
#endif
#if defined(CERNLIB_FQNEEDPK)
      IF (IUPAKI.NE.0)             GO TO 31

C--       Exchange Data Format needing to be unpacked

      LIN = LIN + IQ(KQSP+LBPARI+3)
#endif
#if (!defined(CERNLIB_QREADFULL))&&(defined(CERNLIB_FQNEEDPK))
      IF (NW4IN.EQ.MAXREI)         GO TO 31
#include "fznwmach.inc"
* Ignoring t=pass
#endif

C----              READ ONE PHYSICAL RECORD

   31 IQ(KQSP+LQFI+33) = IQ(KQSP+LQFI+33) + 1
      NBLK = IQ(KQSP+LQFI+22) + 1
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +      WRITE (IQLOG,9031) NBLK,NW4IN,NWMIN,NFASTI,NRSKIP
 9031 FORMAT (1X/' FZIPHR-  Reading Block',I7,
     F', NW32,NWmach,NRfast,NRskip=',4I6)

#endif
#if defined(CERNLIB_FZLIBC)
      IF (IACMOI.EQ.2)  THEN
          NWR = NWMREC
          CALL CFGET (IADOPI, MEDIUI, NWMREC, NWR, LQ(LIN), ISW)
          IF (ISW.EQ.-1)           GO TO 841
          IF (ISW.NE.0)            GO TO 843
          IF (NWR.EQ.NWMREC)  NW4IN = MAXREI
          GO TO 39
        ENDIF

#endif
#if defined(CERNLIB_FZCHANNEL)
      IF (IACMOI.EQ.3)  THEN
          CALL JUMPST (IADOPI)
          ICODE     = 0
          IQUEST(1) = LUNI
          IQUEST(2) = NWMREC
          IQUEST(3) = ISTENI
          IQUEST(4) = 0
          IQUEST(5) = MEDIUI - 4
          IQUEST(6) = NWMREC
          CALL JUMPX2 (LQ(LIN),ICODE)
          ISW = IQUEST(1)
          IF (ISW.LT.0)            GO TO 841
          IF (ISW.NE.0)            GO TO 843
          NWR = IQUEST(2)
          IF (NWR.EQ.NWMREC)  NW4IN = MAXREI
          GO TO 39
        ENDIF

#endif
#if defined(CERNLIB_FZFORTRAN)
#include "fziphr32.inc"
   36 NWR = NWMIN
      CALL XINBF (LUNI,LQ(LIN),NWR)
      IF (NWR.EQ.0)                GO TO 841
      ISW = -NWR
      IF (NWR.LT.0)                GO TO 843
      NWR = MIN (NWR,NWMIN)
#endif
   39 IQ(KQSP+LQFI+22) = NBLK
      NWRDAI = NWRDAI + MAXREI
      NW4USE = NW4IN

C----              UNPACK / BYTE-SWOP

#if defined(CERNLIB_FQNEEDCV)
      IF (IUPAKI.NE.0)             GO TO 47
      IF (NRSKIP.NE.0)   NW4USE = 8
      IF (IFLAGI.EQ.-1)  NW4USE = 8

#endif
#if defined(CERNLIB_FQNEEDCV)
#include "fziphr42.inc"
#endif

C--                Short/full dump of record read

   47 CONTINUE
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  CALL FZIDUM (LQ(L4STAI),NW4USE)
#endif

C-----------------------------------------------------------
C----              FAST RECORD EXPECTED
C-----------------------------------------------------------

      IF (NFASTI.EQ.0)             GO TO 61
      IF (LQ(L4STAI).EQ.MCCW1)     GO TO 54
   52 NFASTI = NFASTI - 1
      IQ(KQSP+LBPARI-5) = NFASTI
      IF (NRSKIP.EQ.0)             GO TO 53

C--       skip record
      N4SKII = N4SKII - MAXREI
      N4RESI = N4RESI - MAXREI
      NRSKIP = NRSKIP - 1
      GO TO 20

C--       deliver record
   53 N4DONI = 0
      N4ENDI = MIN (N4RESI,MAXREI)
      IFLAGI = 0
#include "zebra/qtrace99.inc"
      RETURN

C----              Unexpected steering record

   54 CALL FZICHH (0, LQ(L4STAI),1)
      IF (IQUEST(1).EQ.0)          GO TO 52
#if defined(CERNLIB_FQNEEDCV)
C--       Unpack all words read
      IF (NW4USE.LT.NW4IN) THEN
          NW4USE = NW4IN
          GO TO 42
        ENDIF
#endif

      IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      N4ENDI = NTLRI
      IF (JMODI.EQ.4)              GO TO 809
      IF (NW4IN.EQ.MAXREI)         GO TO 808

C--         recover start/end-of-run in partially read record

      IF (JMODI.EQ.0)              GO TO 806
      IF (JMODI.GT.2)              GO TO 806
      IF (NTLRI+3.GT.NW4IN)        GO TO 805

      LRCUR = L4STAI + NTLRI + 1
      LRLG  = LQ(LRCUR-1)
      LRTP  = LQ(LRCUR)
      IF (LRLG.LT.1)               GO TO 806
      IF (LRLG.GT.401)             GO TO 806
      IF (LRTP.NE.1)               GO TO 806
      IF (LQ(LRCUR+1).LT.0)        GO TO 806
      LRLG  = MIN (LRLG, NW4IN-NTLRI-2)
      NUSED = NTLRI + LRLG + 2
      LQ(LRCUR-1) = LRLG
      LRCUR = LRCUR + LRLG + 1
      LQ(LRCUR)   = MAXREI - NUSED - 1
      LQ(LRCUR+1) = 5
      NFSTI  = 0
      GO TO 807

C-----------------------------------------------------------
C----              STEERING RECORD EXPECTED
C-----------------------------------------------------------

   61 CALL FZICHH (0, LQ(L4STAI),IFLAGI)
      IF (IQUEST(1).NE.0)          GO TO 71
      N4ENDI = NTLRI

      IF (IFLAGI.LT.0)             GO TO 73

      JREX = IQ(KQSP+LBPARI-7)
      IF   (JREX.NE.0)  THEN
          IF (JRECI.NE.0)    THEN
              IF (JRECI.NE.JREX)   GO TO 804
            ENDIF
        ENDIF

      N4DONI = 8
      IF (N4ENDI.EQ.0)  N4ENDI= MAXREI

   62 IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      IF (JRECI.NE.0)   JRECI = JRECI + 1
      IQ(KQSP+LBPARI-7) = JRECI
      IQ(KQSP+LBPARI-5) = NFSTI
      IFLAGI = 0
      GO TO 999

C----              Recover to next steering record

   71 IF (IQUEST(1).EQ.3)          GO TO 802
      IF (IFLAGI.EQ.-1)            GO TO 20
      GO TO 801

C--                Recovery to this steering record

   73 IF (NTLRI.EQ.0)              GO TO 20
#if defined(CERNLIB_FQNEEDCV)
      IF (NW4USE.LT.NW4IN) THEN
          NW4USE = NW4IN
          GO TO 42
        ENDIF
#endif
      N4DONI = NTLRI
      GO TO 62

C-----------------------------------------------------------
C-                 ERROR CONDITIONS
C-----------------------------------------------------------


C-    JERROR = 201  Block header faulty
  801 JERROR = 201
      GO TO 817

C-    JERROR = 202  Block size does not match expectation
  802 JERROR = 202
      IQUEST(14) = MAXREI
      IQUEST(15) = NWRI
      NWERR = 2
      GO TO 817

C-    JERROR = 204  Break in block sequence number
  804 JERROR = 204
      JRETCD = 5
      IQUEST(14) = JREX
      IQUEST(15) = JRECI
      NWERR  = 2
      GO TO 811

C-    JERROR = 205  Fast burst stopped by unusable start/end-of-run
  805 JERROR = -1

C-    JERROR = 206  Fast burst stopped by unusable steering block
  806 JERROR = JERROR - 1
      N4ENDI = 0

C-    JERROR = 207  Fast burst stopped by usable start/end-of-run
C-                             in unusable steering block
  807 JERROR = JERROR - 1

C-    JERROR = 208  Fast burst stopped by usable steering block
  808 JERROR = 208 + JERROR
      IQUEST(14) = NTLRI
      IQUEST(15) = LQ(L4STAI+8)
      IQUEST(16) = LQ(L4STAI+9)
      NWERR  = 3
      JRETCD = 5
      GO TO 811

C-    JERROR = 209  Emergency stop block
  809 JERROR = 209
      JRETCD = 8
      N4ENDI = 0

  811 IQ(KQSP+LBPARI-7) = 0
      IQ(KQSP+LBPARI-6) = 0
      IQ(KQSP+LBPARI-5) = NFSTI
      IQ(KQSP+LBPARI-1) = N4ENDI
      GO TO 819

  817 JRETCD = 6
  818 IQ(KQSP+LBPARI-1) = 0
  819 IQ(KQSP+LBPARI-9)= -1
  820 IFLAGI = 1
      GO TO 999

C--                EoF seen

  841 JRETCD = 1
      IQ(KQSP+LBPARI-7) = 0
      IQ(KQSP+LBPARI-5) = 0
      IQ(KQSP+LBPARI-1) = 0
      GO TO 820

C--                Read error

  843 JRETCD = 7
      JERROR = 215
      NWERR  = 1
      IQUEST(14) = ISW
      GO TO 818
      END
*      ==================================================
#include "zebra/qcardl.inc"
