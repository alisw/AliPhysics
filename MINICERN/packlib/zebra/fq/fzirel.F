*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:43  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIREL

C-    RELOCATE DATA-STRUCTURE READ
C-    USE THE MEMORY OCCUPATION TABLE READY
C-    AND THE RELOCATION VECTOR READ INTO LQ(LQTA+NWTABI)
C-    CALLED FROM FZIN

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    LADESV(6)
#if defined(CERNLIB_QMVDS)
      SAVE         LADESV
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZIR, 4HEL   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIREL /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIREL  ')
#endif
      DATA  LADESV / 6, 5*0 /


#include "zebra/qtrace.inc"

      IF (NWTABI.EQ.0)             GO TO 61

C------            Ready the relocation table

      LPUT  = LQTA
      LTAKE = LQTA + NWTABI

C--                Loop for all segments in the memory occupation tb.

      LMT = LQMTA
   22 IF (LQ(LMT+1).NE.0)          GO TO 24

C--                Skipped segment

      NWSG  = LQ(LMT+3)
   23 IF (NWSG.GE.0)               GO TO 29
      IF (LTAKE.GE.LQTE)           GO TO 731
      NWSG  = NWSG + (LQ(LTAKE+1)-LQ(LTAKE))
      LTAKE = LTAKE + 2
      GO TO 23

C--                Accepted segment

   24 LSTA  = LQ(LMT+3)
      LEND  = LQ(LMT+4)
      NWSG  = LSTA - LEND
      NREL  = 0
      LE    = LSTA
   25 IF (LTAKE.GE.LQTE)           GO TO 731
      LA    = LQ(LTAKE)
      NREL  = NREL - (LA-LE)
      LE    = LQ(LTAKE+1)

      LQ(LPUT)   = LA
      LQ(LPUT+1) = LE
      LQ(LPUT+2) = NREL
      LQ(LPUT+3) = 0
      LTAKE = LTAKE + 2
      LPUT  = LPUT  + 4
      NWSG  = NWSG + (LE-LA)
      IF (NWSG.LT.0)               GO TO 25
   29 IF (NWSG.NE.0)               GO TO 732
      LMT = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 22
      IF (LTAKE.NE.LQTE)           GO TO 733
      LQTE = LPUT

      LQ(LQTE)   =  LQ(LQTE-3)
      LQ(LQTA-1) =  LQ(LQTA)
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.4)
     +  WRITE (IQLOG,9167) LENTRI,(LQ(J),J=LQTA,LQTE-1)
 9167 FORMAT (' FZIREL-  Relocation Table, LENTRY before=',I10/
     F (15X,3I9,I4))
#endif

C----              Relocate the bank links

      IQFLIO = 7
      CALL MZRELB
      IF (IQFLIO.LT.0)             GO TO 734

C--                Relocate the entry link

      LADESV(2) = LOCF(LENTRI) - LQSTOR
      LADESV(3) = LADESV(2) + 1
      LADESV(5) = IQLETT(9)
      LADESV(6) = IQLETT(15)
      CALL MZRELL (LADESV)
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.4)  WRITE (IQLOG,9037) LENTRI
 9037 FORMAT (10X,'LENTRY after=',I10)
#endif
      LQ(KQS+LENTRI+1) = 0
      LQ(KQS+LENTRI+2) = 0
      GO TO 999

C------            Chain banks into one linear structure

   61 CALL FZILIN
      IF (IQFOUL.NE.0)             GO TO 734
      LENTRI = IQUEST(1)
#include "zebra/qtrace99.inc"
      RETURN

C-------------------------------------------------
C-                 ERROR CONDITIONS
C-------------------------------------------------

C----              BAD DATA

C-    JERROR = 34  bank chaining clobbered in the input data
  734 JERROR = 34
      IQUEST(14)= IQLN
      NWERR  = 1
      GO TO 739

C-    JERROR = 33  ends of segment and rel. tables do not match
  733 JERROR = 33
      IQUEST(14)= LTAKE
      IQUEST(15)= LQTE
      NWERR  = 2
      GO TO 739

C-    JERROR = 32  segment limit does not match a rel. table entry
  732 JERROR = 32
      IQUEST(14)= NWSG
      NWERR  = 1
      GO TO 739

C-    JERROR = 31  segment table tries to overshoot rel. table
  731 JERROR = 31
      IQUEST(14)= NWSG
      NWERR  = 1
  739 JRETCD = 5
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
