*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:38  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZDACC)
      SUBROUTINE FZIPHD

C-    Read next physical record in direct-access mode
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
#include "fzdaeof.inc"
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIP, 4HHD   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIPHD /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIPHD  ')
#endif

*      Declaratives, DATA
#include "fziphrd2.inc"


#include "zebra/qtrace.inc"

      NWMREC = IQ(KQSP+LBPARI+1)
      NRSKIP = 0

      JRECGO = IQ(KQSP+LQFI+33)
      JRECLL = JRECGO
      IF (JRECGO.GE.0)             GO TO 22

C--                Random access

      JRECNO = IQ(KQSP+LQFI+31)
      NFASTI = 0
      GO TO 27

C--                Sequential access

   22 NFASTI = IQ(KQSP+LBPARI-5)
      JRECNO = JRECGO + 1
      N4SKIP = MIN (N4SKII,N4RESI)
      IF (JRECGO.EQ.0)  JRECLL = -7

      IF (N4SKIP.GE.MAXREI)  THEN
        IF (NFASTI.GT.0)  THEN
            NRSKIP = MIN (NFASTI, N4SKIP / MAXREI)
            JRECNO = JRECNO + NRSKIP
            NFASTI = NFASTI - NRSKIP
          ENDIF
        ENDIF

      IQ(KQSP+LBPARI-6) = NFASTI

   27 LIN = L4STAI
#if defined(CERNLIB_FQNEEDPK)
      IF (IUPAKI.NE.0)             GO TO 31

C--       Exchange Data Format needing to be unpacked

      LIN = LIN + IQ(KQSP+LBPARI+3)
#endif

C----              READ ONE PHYSICAL RECORD

   31 NBLK = IQ(KQSP+LQFI+22) + 1
      IQ(KQSP+LQFI+33) = JRECNO
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +      WRITE (IQLOG,9031) JRECNO,MAXREI,NWMREC,NFASTI,NRSKIP
 9031 FORMAT (1X/' FZIPHD-  Reading Block #',I7,
     F', NW32,NWmach,NRfast,NRskip=',4I6)
#endif
#if defined(CERNLIB_FZDACCL)
      IF (IACMOI.EQ.2)  THEN
          IF (JRECNO.NE.JRECLL+1)  THEN
              CALL CFSEEK (IADOPI, MEDIUI, NWMREC, JRECNO-1, ISW)
              IF (ISW.NE.0)        GO TO 847
            ENDIF

          NWR = NWMREC
          CALL CFGET (IADOPI, MEDIUI, NWMREC, NWR, LQ(LIN), ISW)
          IF (ISW.NE.0)            GO TO 842
          JRECLL = JRECNO
          GO TO 39
        ENDIF

#endif
#if defined(CERNLIB_FZDACCH)
      IF (IACMOI.EQ.3)  THEN
          CALL JUMPST (IADOPI)
          ICODE     = 0
          IQUEST(1) = LUNI
          IQUEST(2) = NWMREC
          IQUEST(3) = ISTENI
          IQUEST(4) = JRECNO
          IQUEST(5) = MEDIUI - 4
          IQUEST(6) = NWMREC
          CALL JUMPX2 (LQ(LIN),ICODE)
          ISW = IQUEST(1)
          IF (ISW.LT.0)            GO TO 841
          IF (ISW.NE.0)            GO TO 843
          NWR = IQUEST(2)
          GO TO 39
        ENDIF

#endif
#if defined(CERNLIB_FZDACCF)
#include "fziphd37.inc"
#endif
   39 IQ(KQSP+LQFI+22) = NBLK
      NWRDAI = NWRDAI + MAXREI
      NW4USE = MAXREI

C----              UNPACK / BYTE-SWOP

#if defined(CERNLIB_FQNEEDCV)
      IF (IUPAKI.NE.0)             GO TO 47

#endif
#if defined(CERNLIB_FQNEEDCV)
#include "fziphd42.inc"
#endif

C--                Short/full dump of record read

   47 CONTINUE
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)  CALL FZIDUM (LQ(L4STAI),NW4USE)
#endif

C-----------------------------------------------------------
C----              FAST RECORD EXPECTED
C-----------------------------------------------------------

      IF (NFASTI.LT.0)             GO TO 84
      IF (NFASTI.EQ.0)             GO TO 61
      IF (LQ(L4STAI).EQ.MCCW1)     GO TO 57
   52 NFASTI = NFASTI - 1
      IQ(KQSP+LBPARI-5) = NFASTI
      IF (NRSKIP.NE.0)  THEN
C--       skip records
          N4SKII = N4SKII - MAXREI*NRSKIP
          N4RESI = N4RESI - MAXREI*NRSKIP
        ENDIF

C--       deliver record
      N4DONI = 0
      N4ENDI = MIN (N4RESI,MAXREI)

      IFLAGI = 0
#include "zebra/qtrace99.inc"
      RETURN

C----              Unexpected steering record

   57 CALL FZICHH (0, LQ(L4STAI),1)
      IF (IQUEST(1).EQ.0)          GO TO 52
      GO TO 81

C-----------------------------------------------------------
C----              STEERING RECORD EXPECTED
C-----------------------------------------------------------

   61 CALL FZICHH (0, LQ(L4STAI),IFLAGI)
      IF (IQUEST(1).NE.0)          GO TO 71
      IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      N4DONI = 8

C--       random access
      IF (JRECGO.LT.0)  THEN
          N4ENDI = IQ(KQSP+LQFI+32)
          IF (N4ENDI.EQ.0)  N4ENDI = NTLRI
          IF (JRECI.NE.0)   JRECI  = JRECI + 1
          GO TO 67
        ENDIF

C--                Sequental access

      N4ENDI = NTLRI
      IF (N4ENDI.EQ.0)  N4ENDI= MAXREI

      IF (IFLAGI.EQ.-1)            GO TO 73

C--       skip records
      IF (NRSKIP.NE.0)  THEN
          N4SKIP = N4SKIP - MAXREI*NRSKIP
          N4RESI = N4RESI - MAXREI*NRSKIP
        ENDIF

      JREX = IQ(KQSP+LBPARI-7)
      IF   (JRECI.NE.0)  THEN
          IF (JREX.NE.0)    THEN
              IF (JREX.NE.JRECI)   GO TO 804
            ENDIF
          JRECI = JRECI + 1
        ENDIF

   67 IQ(KQSP+LBPARI-7) = JRECI
      IQ(KQSP+LBPARI-5) = NFSTI
      IFLAGI = 0
      GO TO 999

C----              Recover to next steering record

   71 IF (IQUEST(1).EQ.3)          GO TO 802
      IF (IQUEST(1).EQ.4)          GO TO 841
      IF (IFLAGI.NE.-1)            GO TO 82
   72 JRECNO = JRECNO + 1
      GO TO 27

C--                Recovery to this steering record

   73 IF (NTLRI.EQ.0)              GO TO 72
      N4DONI = NTLRI
      GO TO 67

C-----------------------------------------------------------
C----              Check skipping passed over a steering record
C-----------------------------------------------------------

   81 IF (NRSKIP.EQ.0)             GO TO 87
      GO TO 83

   82 IF (NRSKIP.EQ.0)             GO TO 801
   83 JRECE  = JRECNO
      JRECNO = JRECGO + 1
      NFASTI = -1
      GO TO 27

   84 CALL FZICHH (0, LQ(L4STAI),-1)
      IF (IQUEST(1).EQ.0)          GO TO 87
      IF (JRECNO.EQ.JRECE)         GO TO 801
      JRECNO = JRECNO + 1
      GO TO 27

   87 N4ENDI = NTLRI
      IF (JMODI.EQ.4)              GO TO 809
      IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      GO TO 808

C-----------------------------------------------------------
C----              ERROR CONDITIONS
C-----------------------------------------------------------

C-    JERROR = 501  Block header faulty
  801 JERROR = 501
      GO TO 817

C-    JERROR = 502  Block size does not match expectation
  802 JERROR = 502
      IQUEST(14) = MAXREI
      IQUEST(15) = NWRI
      NWERR = 2
      GO TO 817

C-    JERROR = 504  Break in block sequence number
  804 JERROR = 504
      JRETCD = 5
      IQUEST(14) = JREX
      IQUEST(15) = JRECI
      NWERR  = 2
      GO TO 811

C-    JERROR = 508  Fast burst stopped by usable steering block
  808 JERROR = 508
      IQUEST(14) = NTLRI
      IQUEST(15) = LQ(L4STAI+8)
      IQUEST(16) = LQ(L4STAI+9)
      NWERR  = 3
      JRETCD = 5
      GO TO 811

C-    JERROR = 509  Emergency stop block
  809 JERROR = 509
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

  840 IF (ISW.EQ.-1)               GO TO 841
*      Extra EOF test
#include "fzdaeofx.inc"
* Ignoring t=pass
      IF (ISW.NE.IOSEOF)           GO TO 843
  841 JRETCD = 1
      IQ(KQSP+LBPARI-7) = 0
      IQ(KQSP+LBPARI-5) = 0
      IQ(KQSP+LBPARI-1) = 0
      GO TO 820

C--                Read error

#if defined(CERNLIB_FZDACCL)
  842 IF (ISW.EQ.-1)               GO TO 841
#endif
  843 JRETCD = 7
C-    JERROR = 515  Read error
      JERROR = 515
      NWERR  = 1
      IQUEST(14) = ISW
      GO TO 818

C--                Seek error

#if defined(CERNLIB_FZDACCL)
  847 JRETCD = 7
C-    JERROR = 514  Seek error
      JERROR = 514
      NWERR  = 1
      IQUEST(14) = ISW
      GO TO 818
#endif
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
