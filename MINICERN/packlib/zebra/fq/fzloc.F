*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:44  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZLOC (LUNP,MODEP)

C-    Locate FZ control-bank for unit LUN
C-    Unpack control information

C-    MODE =  -2 or +2 : connect as output unit
C-            -1 or +1 : connect as input  unit
C-                   0 : no connection (for FZFILE, FZHOOK, FZLIMI, etc)
C-         if > 0 : unit is required for the given I/O mode
C-         if < 0 : wanting to 'end' the unit (for FZENDx)
C-                  if IQUEST(1) < 0 : connect only if active
C-                                     (loop for all units)
C-                               = 0 : connect if possible
C-                                     (for the particular unit)

#include "zebra/zstate.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcf.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9), MODEP(9)
      DIMENSION    PILX(4)
      EQUIVALENCE (PILX(1),IPILX(1))
#if defined(CERNLIB_QMVDS)
      SAVE         CHDATA
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZLO, 4HC    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZLOC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZLOC   ')
#endif
      DATA  CHDATA / 12345.0 /

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      LUN  = LUNP(1)
      MODE = MODEP(1)
      MODA = IABS(MODE)
      MODS = IQUEST(1)

      LQFF = LQFS
   14 IF (LQFF.NE.0)  THEN
          IF (IQ(KQSP+LQFF-5).EQ.LUN)  GO TO 21
          LQFF = LQ(KQSP+LQFF)
          GO TO 14
        ELSE
          LUNF  = 0
          IACTVF= 0
          IF (MODE.GT.0)           GO TO 94
          RETURN
        ENDIF

C--                Bank for LUN found

   21 MSTATF = IQ(KQSP+LQFF)
      MEDIUF = JBYT (MSTATF,1,3)
      IFIFOF = JBYT (MSTATF,4,3)
      IDAFOF = JBIT (MSTATF,7)
      IACMOF = JBYT (MSTATF,8,3)
      IUPAKF = JBIT (MSTATF,16)
      IADOPF = IQ(KQSP+LQFF+1)
      IACTVF = IQ(KQSP+LQFF+2)
      INCBPF = IQ(KQSP+LQFF+3)
      LOGLVF = IQ(KQSP+LQFF+4)
      MAXREF = IQ(KQSP+LQFF+5)
      LUNF   = LUN
      IF   (MODA-1)          79, 31, 51

C--                Connect as input unit

   31 IF (IACTVF.LT.1)             GO TO 34
      IF (IACTVF.GE.8)             GO TO 34
   32 MSTATI = MSTATF
      MEDIUI = MEDIUF
      IFIFOI = IFIFOF
      IDAFOI = IDAFOF
      IACMOI = IACMOF
      IUPAKI = IUPAKF
      IADOPI = IADOPF
      IACTVI = IACTVF
      INCBPI = INCBPF
      LOGLVI = LOGLVF
      MAXREI = MAXREF
      LUNI   = LUNF
      LQFI   = LQFF
      GO TO 79

C--   connection required
   34 IF (MODE.LT.0)               GO TO 36
      IF (JBIT(MSTATF,11).EQ.0)    GO TO 93
      IF (IACTVF.EQ.0)             GO TO 37
      IF (IACTVF.EQ.8)             GO TO 37
      IF (IACTVF.NE.18)            GO TO 92
      IACTVF = 0
      GO TO 37

C--   connect if possible
   36 IF (MODS.LT.0)               GO TO 79
      IF (JBIT(MSTATF,11).EQ.0)    GO TO 79
      IF (IACTVF.EQ.0)             GO TO 37
      IF (IACTVF.EQ.8)             GO TO 37
      IF (IACTVF.NE.18)            GO TO 79

C--   clear buffer parameters
   37 IF (IFIFOF.EQ.0)             GO TO 32
      CALL VZERO (IQ(KQSP+LQFF+40),INCBPF-40)
      GO TO 32


C--                Connect as output unit

   51 IF (IACTVF.LT.11)            GO TO 54
      IF (IACTVF.EQ.18)            GO TO 54
   52 MSTATX = MSTATF
      MEDIUX = MEDIUF
      IFIFOX = IFIFOF
      IDAFOX = IDAFOF
      IACMOX = IACMOF
      IUPAKX = IUPAKF
      IADOPX = IADOPF
      IACTVX = IACTVF
      INCBPX = INCBPF
      LOGLVX = LOGLVF
      MAXREX = MAXREF
      LUNX   = LUNF
      LQFX   = LQFF
      GO TO 79

C--   connection required
   54 IF (MODE.LT.0)               GO TO 56
      IF (JBIT(MSTATF,12).EQ.0)    GO TO 93
      IF (IACTVF.EQ.10)            GO TO 57
      IF (IACTVF.NE.0)             GO TO 91
      GO TO 57

C--   connect if possible
   56 IF (MODS.LT.0)               GO TO 79
      IF (JBIT(MSTATF,12).EQ.0)    GO TO 79
      IF (IACTVF.EQ.10)            GO TO 57
      IF (IACTVF.NE.0)             GO TO 79

C--   clear buffer parameters
   57 IF (IFIFOF.NE.0)  CALL VZERO (IQ(KQSP+LQFF+40),INCBPF-40)
       PILX(1) = CHDATA
      IPILX(2) = 10000.0 * QVERSN + .2
      IPILX(3) = 0
      IPILX(4) = 0
      GO TO 52

   79 RETURN

C-------------------------------------------------
C-                 trouble
C-------------------------------------------------

C--                Write after read without FZEND

   91 NQCASE = -1

C--                Read after write without REWIND

   92 NQCASE = NQCASE + 2
      NQFATA = 3
      IQUEST(13) = IACTVF
      GO TO 90

C--                Permission fault

   93 NQCASE = 3
      NQFATA = 3
      IQUEST(13) = JBYT (MSTATF,11,2)
      GO TO 90

C--                File not opened

   94 NQCASE = 4
      NQFATA = 2
   90 IQUEST(11) = LUN
      IQUEST(12) = MODE
#include "zebra/qtrace.inc"
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
