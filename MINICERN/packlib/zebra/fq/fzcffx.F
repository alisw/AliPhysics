*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZCFFX

C-    Copy table + bank material for input file format X,
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
      DATA  NAMESR / 4HFZCF, 4HFX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZCFFX  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZCFFX  ')
#endif

#include "zebra/qtraceq.inc"


      NWDO   = NWBKX
      ISTTAB = IQUEST(7)
#if defined(CERNLIB_FZFFNAT)
      IF (IFIFOX.EQ.0)             GO TO 41
#endif

C-------------------------------------------------
C-                 Output file format exchange
C-------------------------------------------------

      IF (ISTTAB.LT.0)  NWDO = NWDO + NWTABX

   24 CALL FZIACN (NWDO,LDATA,NWACC)
      IF (IFLAGI.NE.0)             GO TO 999

      CALL FZOTRN (LQ(LDATA),NWACC)
      NWDO = NWDO - NWACC
      IF (NWDO.NE.0)               GO TO 24

   29 CALL FZIREC
      IF (N4RESI.NE.0)             GO TO 991

      IF (IFIFOX.NE.0)  CALL FZOREC

#include "zebra/qtrace99.inc"
      RETURN

C-------------------------------------------------
C-                 Output file format native
C-------------------------------------------------

#if defined(CERNLIB_FZFFNAT)
   41 MINREC = (4*MAXREX) / 5
      LAST   = 0
      IF (ISTTAB.GE.0)             GO TO 61

C----              Copy the table

      IDX(2) = 4
      NWDO = NWTABX

   44 NWR  = MIN (NWDO, MAXREX)
      CALL FZIACN (NWR,LDATA,NWACC)
      IF (IFLAGI.NE.0)             GO TO 999
      NWDO = NWDO - NWACC
      IF (NWDO.EQ.0)  THEN
          IF (LAST.NE.0) IDX(2)=8
          CALL FZON1 (LQ(LDATA),NWACC)
          GO TO 58
        ENDIF

      IF (NWACC.GE.MINREC)  THEN
          CALL FZON1 (LQ(LDATA),NWACC)
          GO TO 44
        ENDIF

      NIN = 0
   46 CALL UCOPY (LQ(LDATA),IQWKTB(NIN+1),NWACC)
      NIN = NIN + NWACC

      NWR  = MIN (NWDO, MAXREX-NIN)
      CALL FZIACN (NWR,LDATA,NWACC)
      IF (IFLAGI.NE.0)             GO TO 999
      NWDO = NWDO - NWACC
      IF (NWDO.EQ.0)  THEN
          IF (LAST.NE.0)  IDX(2)=8
        ELSE
          IF (NIN+NWACC.LT.MINREC) GO TO 46
        ENDIF

      CALL FZON2 (IQWKTB,NIN,LQ(LDATA),NWACC)
      IF (NWDO.NE.0)               GO TO 44

   58 IF (LAST.NE.0)               GO TO 999
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
#endif
C----              Error handling

C-    JERROR = 455  bank material does not end exactly with LR
  991 JERROR = 455
      JRETCD = 5
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
