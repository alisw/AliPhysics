*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:30  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZDACC)
      SUBROUTINE FZIDAT (LUNP,IXDIVP,LSUPP,JBIASP)

C-    Find and read the Direct Access Table

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9),IXDIVP(9),LSUPP(9),JBIASP(9)
      DIMENSION    MDSADR(4)
      CHARACTER    CHOPT*2
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZID, 4HAT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIDAT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIDAT  ')
#endif

#include "zebra/q_jbit.inc"

#include "zebra/qtrace.inc"

      LUNNI = LUNP(1)
      CALL FZLOC (LUNNI,1)
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.2)  WRITE (IQLOG,9002) LUNI
 9002 FORMAT (' FZIDAT-  called for LUN=',I4)
#endif

      IFLUPD = 7
      IFLNOR = 0
      LBPARI = LQFI + INCBPI
      NWPREC = IQ(KQSP+LBPARI+1)

C----              Read the first logical record

      CHOPT = 'SF'
      ITER  = 0
      MDSADR(1) = 1
      MDSADR(2) = 0

   24 CALL FZINXT (LUNI,MDSADR(1),MDSADR(2))

   26 NUH = 2
      CALL FZIN (LUNI, IXDIVP, LSUPP, JBIASP, CHOPT, NUH, MDSADR)
      MDSADR(NUH+1) = 0
      MDSADR(NUH+2) = 0

C--       ITER = 0  read first LR of the file, is it DAT forward ref.?
C-               1  read the last LR of the file, is it Zebra EoF ?
C-               2  hunt for DAT record
C-               3  select DAT record
C-               4  read the data of the DAT record

      IF (ITER.EQ.4)               GO TO 67
      IF (ITER.EQ.3)               GO TO 64
      IF (ITER.EQ.2)               GO TO 54
      IF (ITER.EQ.1)               GO TO 44

C--       look at the first record
      IF (IQUEST(1).GE.4)          GO TO 802
      IF (IQUEST(1).NE.0)          GO TO 34
      IF (IPILI(3) .NE.2)          GO TO 34
      IF (MDSADR(1).NE.0)          GO TO 61

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9032) LUNI
 9032 FORMAT (' FZIDAT.  LUN=',I4,' DaT forward ref. record not filled')
#endif
      GO TO 41

   34 IFLUPD = 0
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9034) LUNI
 9034 FORMAT (' FZIDAT.  LUN=',I4,
     F' does not start with DaT forward ref. record')
#endif

C------            Get the last logical record

   41 ITER  = 1
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9041)
 9041 FORMAT (19X,'try direct read of last record')
#endif
#if defined(CERNLIB_FZDACCL)
      IF (IACMOI.EQ.2)  THEN
          CALL CFSIZE (IADOPI, MEDIUI,NWPREC,JRECLL,ISTAT)
          GO TO 42
        ENDIF
#endif
      CALL FFSIZE (LUNI,NWPREC,JRECLL,ISTAT)

   42 IF (ISTAT.NE.0)  JRECLL = 0
      IF (JRECLL.LE.0)             GO TO 51
      MDSADR(1) = JRECLL
      MDSADR(2) = 0
      GO TO 24

C--       look at the last logical record
   44 IF (IQUEST(1).GE.4)          GO TO 51
      IF (IQUEST(1).LT.0)          GO TO 51
      IF (IQUEST(1).NE.3)          GO TO 26
      IF (MDSADR(1).NE.0)          GO TO 62

#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9046) LUNI
 9046 FORMAT (' FZIDAT.  LUN=',I4,' Zebra EoF does not point to DaT')
#endif
      IF (MDSADR(2).EQ.0)          GO TO 801

C------            Hunt for the DAT bank or Zebra EoF

   51 ITER   = 2
      CHOPT  = 'SH'
      JRECLL = MAX (JRECLL-25,1)
      LDSADR = -1
      NFTOL  = 4
      IF (JRECLL.GE.11)  NFTOL  = 12
      MDSADR(1) = JRECLL
      MDSADR(2) = 0
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9051) LUNI,JRECLL
 9051 FORMAT (' FZIDAT.  LUN=',I4,' cannot get DaT adr from Zebra EoF'
     F/19X,'hunt for it starting at record',I9)
#endif
      GO TO 24

   54 IF (IQUEST(1).EQ.0)          GO TO 65
      IF (IQUEST(1).EQ.3)          GO TO 55
      IF (IQUEST(1).GE.4)          GO TO 57
      IF (IQUEST(1).LT.0)  THEN
          NFTOL = NFTOL - 1
          IF (NFTOL.LT.0)          GO TO 801
        ELSE
          LDSADR = -1
        ENDIF
      GO TO 26

   55 IF (MDSADR(1).NE.0)          GO TO 62
      LDSADR = MDSADR(2)
      GO TO 26

   57 IF (JRECLL.LE.1)             GO TO 801
      IF (LDSADR.EQ.0)             GO TO 801
      JRECLL = 0
      GO TO 51

C------            Read the DAT bank

   61 IFLUPD = 0
      IFLNOR = -2
   62 ITER   = 3
      GO TO 24

   64 IF (IQUEST(1).NE.0)          GO TO 801
      IF (IPILI(3).NE.1)           GO TO 801
   65 IQ(KQSP+LQFI+34) = IQUEST(5)
      IQ(KQSP+LQFI+35) = IQUEST(6)
      CHOPT = 'A '
      ITER  = 4
      GO TO 26

   67 IF (IQUEST(1).NE.0)          GO TO 801
      IF (IQUEST(13).EQ.0)         GO TO 801

C--       update the DAT forward reference record
      IF (IFLUPD.EQ.0)             GO TO 71
      IF (JBIT(MSTATI,12).EQ.0)    GO TO 71
      IQ(KQSP+LQFI+2) = 0
      CALL FZUDAT (LUNP,1)

   71 IQ(KQSP+LQFI+2) = 2
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI+IFLNOR.GE.0)  WRITE (IQLOG,9072) LUNI
 9072 FORMAT (' FZIDAT.  LUN=',I4,' rewind')
#endif
      CALL FZENDI (LUNP,'IQ')
#include "zebra/qtrace99.inc"
      RETURN

C----              Errors

  801 IQUEST(1) = -1
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.-2)  WRITE (IQLOG,9098) LUNI
 9098 FORMAT (' FZIDAT.  LUN=',I4,' DaT not found !!!')
#endif
      GO TO 71

  802 IQUEST(1) = -2
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVI.GE.-2)  WRITE (IQLOG,9099) LUNI
 9099 FORMAT (' FZIDAT.  LUN=',I4,' file is empty !!!')
#endif
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
