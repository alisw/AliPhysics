*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:38  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZIPHA

C-    Read next physical record image in ALFA exchange mode
C-    Service routine to FZIN, called only via FZIREC

C-    Input :      IFLAGI = 0  normal read
C-                         -1  recover to next steering block

C-                 N4SKII      is used for rapid skip of fast blocks

C-    Output :     IFLAGI = 0  all is well
C-                             otherwise : ready for re-start

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
#include "fzhci.inc"
C--------------    End CDE                             --------------

*      Declaratives, DIMENSION etc.
#include "fzstamp.inc"

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIP, 4HHA   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIPHA /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIPHA  ')
#endif


#include "zebra/qtrace.inc"

C----              Skip to next steering block

      IF (IFLAGI.LT.0)  THEN
          NFASTI = 0
          NRSKIP = 1
          GO TO 20
        ENDIF

C----              Dedide skip of fast records

      N4SKIP = MIN (N4SKII,N4RESI)
      NFASTI = IQ(KQSP+LBPARI-5)
      NRSKIP = 0
      IF (N4SKIP.LT.MAXREI)        GO TO 20
      IF (NFASTI.EQ.0)             GO TO 20
      NRSKIP = MIN (NFASTI, N4SKIP / MAXREI)

   20 NW4IN  = MAXREI
      IQ(KQSP+LBPARI-6) = NFASTI


C----              Read one physical record

   31 NBLK = IQ(KQSP+LQFI+22) + 1
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.3)
     +      WRITE (IQLOG,9031) NBLK,NW4IN,NFASTI,NRSKIP
 9031 FORMAT (1X/' FZIPHA-  Reading Block',I7,
     F', NW32,NRfast,NRskip=',4I6)
#endif

      CALL FZIASC (NRSKIP)
      IF    (IQUEST(1))      841, 33, 818
   33 JFASTR = IQUEST(92)

      IQ(KQSP+LQFI+22) = NBLK
      NWRDAI = NWRDAI + MAXREI

      IF (JFASTR.EQ.0)  THEN
          LQ(L4STAI)   = MCCW1
          LQ(L4STAI+1) = MCCW2
          LQ(L4STAI+2) = MCCW3
          LQ(L4STAI+3) = MCCW4
        ENDIF

C--                Short/full dump of record read

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.LT.3)             GO TO 51
      IF (JFASTR.NE.0)  THEN
          IF (NRSKIP.NE.0)         GO TO 51
        ENDIF

      CALL FZIDUM (LQ(L4STAI),NW4IN)
#endif

C-----------------------------------------------------------
C----              Fast record expected
C-----------------------------------------------------------

   51 IF (NFASTI.EQ.0)             GO TO 61
      IF (JFASTR.EQ.0)             GO TO 54
      NFASTI = NFASTI - 1
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

C--                Unexpected steering record

   54 CALL FZICHH (0, LQ(L4STAI),0)
      IF (IQUEST(1).NE.0)          GO TO 807

      N4ENDI = NTLRI
      IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      GO TO 808

C-----------------------------------------------------------
C----              Steering record expected
C-----------------------------------------------------------

   61 IF (JFASTR.NE.0)             GO TO 71
      CALL FZICHH (0, LQ(L4STAI),IFLAGI)
      IF (IQUEST(1).NE.0)          GO TO 72
      N4ENDI = NTLRI

      IF (IFLAGI.LT.0)             GO TO 73

      N4DONI = 8
      IF (N4ENDI.EQ.0)  N4ENDI=MAXREI

   62 IQ(KQSP+LQFI+23) = IQ(KQSP+LQFI+23) + 1
      IF (JRECI.NE.0)  JRECI = JRECI + 1
      IQ(KQSP+LBPARI-7) = JRECI
      IQ(KQSP+LBPARI-5) = NFSTI
      IFLAGI = 0
      GO TO 999

C--                Recover to next steering record

   71 IF (IFLAGI.EQ.-1)            GO TO 20
      GO TO 803

   72 IF (IQUEST(1).EQ.3)          GO TO 802
      IF (IFLAGI.EQ.-1)            GO TO 20
      GO TO 801

C--                Recovery to this steering record

   73 IF (NTLRI.EQ.0)              GO TO 20
      N4DONI = NTLRI
      GO TO 62

C-----------------------------------------------------------
C-                 ERROR CONDITIONS
C-----------------------------------------------------------


C-    JERROR = 301  Block header faulty
  801 JERROR = 301
      GO TO 817

C-    JERROR = 302  Block size does not match expectation
  802 JERROR = 302
      IQUEST(14) = MAXREI
      IQUEST(15) = NWRI
      NWERR = 2
      GO TO 817

C-    JERROR = 303  Unexpected fast record
  803 JERROR = 303
      GO TO 817

C-    JERROR = 307  Unexpected and faulty steering block
  807 JERROR = 307
      GO TO 817

C-    JERROR = 308  Unexpected but valid steering block
  808 JERROR = 308
      IQUEST(14) = NTLRI
      IQUEST(15) = LQ(L4STAI+8)
      IQUEST(16) = LQ(L4STAI+9)
      NWERR  = 3
      JRETCD = 5

      IQ(KQSP+LBPARI-7) = 0
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
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
