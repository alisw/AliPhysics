*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:15  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZFFNAT)
      SUBROUTINE FZCFFN

C-    Copy table + bank material for input file format native,
C-    subsidiary to FZCOPY

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzct.inc"
#include "zebra/mzcwk.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcseg.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LRTYP,IDI(2))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZCF, 4HFN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZCFFN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZCFFN  ')
#endif

#include "zebra/qtrace.inc"


      NWDO   = NWBKX
      ISTTAB = IQUEST(7)
      IF (IFIFOX.EQ.0)             GO TO 41

C-------------------------------------------------
C-                 Output file format exchange
C-------------------------------------------------

C----              Copy early table words

      IF    (ISTTAB)               23, 21, 24
   21 CALL FZOTRN (LQ(LQTA), NWTABX)
      GO TO 24

C----              Copy table and bank material

   23 NWDO = NWTABX + NWBKX
   24 NEXT = -7
      LAST = 7
      JIN  = 1

   31 NWR  = NQWKTT
      CALL XINBS (LUNI,LRTYP,1,IQWKTB(JIN),NWR)
      IF (NWR.EQ.0)                GO TO 411
      IF (NWR.LT.0)                GO TO 751
      NRECAI   = NRECAI + 1
      NWRDAI   = NWRDAI + NWR + 3
      IF (LRTYP.EQ.1)              GO TO 424
      IF (LRTYP.LT.4)              GO TO 733
      IF (LRTYP.EQ.5)              GO TO 31
      IF (LRTYP.EQ.6)              GO TO 31
      IF (LRTYP.GE.9)              GO TO 733

      NWDO = NWDO - NWR
      IF (NWDO.LT.0)               GO TO 734
      IF (NWDO.NE.0)  THEN
          IF (LRTYP.EQ.8)          GO TO 735
        ELSE
          IF (LAST.NE.0)  THEN
              IF (LRTYP.NE.8)      GO TO 736
            ENDIF
        ENDIF
      IF (NEXT.EQ.0)               GO TO 46

   36 CALL FZOTRN (IQWKTB,NWR)
      IF (NWDO.NE.0)               GO TO 31

      CALL FZOREC
#include "zebra/qtrace99.inc"
      RETURN

C-------------------------------------------------
C-                 Output file format native
C-------------------------------------------------

C----              Copy early table words

   41 IDX(2) = 4
      MINREC = (4*MAXREX) / 5
      LAST = 0
      NEXT = 0
      IF    (ISTTAB)               43, 42, 61
   42 CALL FZON1 (LQ(LQTA), NWTABX)
      GO TO 61

C----              Copy the table material

   43 NWDO = NWTABX

   44 JIN  = 1
      GO TO 31

   46 JIN = JIN  + NWR
      IF (NWDO.NE.0)  THEN
          IF (JIN.LT.NQWKTT)       GO TO 31
        ENDIF

      JEX = 1
   47 NEX = MIN (JIN-JEX, MAXREX)
      IF (NWDO.EQ.0)  THEN
          IF (LAST.NE.0)  THEN
              IF (JEX+NEX.EQ.JIN)  IDX(2) = 8
            ENDIF
          GO TO 48
        ENDIF

      IF (NEX.GE.MINREC)           GO TO 48
      CALL UCOPY (IQWKTB(JEX),IQWKTB,NEX)
      JIN = NEX + 1
      GO TO 31

   48 CALL FZON1 (IQWKTB(JEX),NEX)
      JEX = JEX + NEX
      IF (JEX.LT.JIN)              GO TO 47
      IF (NWDO.NE.0)               GO TO 44

      IF (LAST.NE.0)               GO TO 999
      IF (IDX(2).NE.4)             GO TO 63

C----              Copy the bank material

   61 IDX(2) = 7
      IF (NQSEG.EQ.0)  THEN
          LAST = 7
          NWDO = NWBKX
          GO TO 44
        ENDIF

      JSEG = 0
   63 JSEG = JSEG + 1
      NWDO = IQSEGD(JSEG)
      IF (JSEG.EQ.NQSEG)  LAST=7
      GO TO 44

C-----------------------------------------------------------
C--                Input errors
C-----------------------------------------------------------

C--                Unexpected end-of-file

  411 JRETCD = 1
      GO TO 999

C--                Unexpected start/end of run

  424 JRETCD = 2
      NWRU   = MIN (78, NWR)
      LFIIOC = LQFI + JAUIOC
      IQ(KQSP+LFIIOC) = NWRU
      CALL UCOPY (IQWKTB,IQ(KQSP+LFIIOC+1),NWRU)
      GO TO 999

C--                BAD DATA

C-    JERROR = 437  emergency stop record seen
  737 JERROR = 437
      JRETCD = 8
      GO TO 999

C-    JERROR = 436  last bank material record needed is not type 8
  736 JERROR = 1

C-    JERROR = 435  premature LR type 8
  735 JERROR = JERROR + 1

C-    JERROR = 434  table or d/s data longer than expected
  734 JERROR = JERROR + 1

C-    JERROR = 433  record of unexpected record type read
  733 JERROR = 433 + JERROR
      IF (LRTYP.EQ.9)              GO TO 737
      JRETCD = 5
      GO TO 999

C--                READ ERROR

  751 JERROR = 451
      NWERR  = 1
      IQUEST(14) = -NWR
      JRETCD = 7
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
