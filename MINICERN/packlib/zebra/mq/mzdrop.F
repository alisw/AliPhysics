*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:24  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:18  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZDROP (IXSTOR,LHEADP,CHOPT)

C-    Drop d/s supported by bank at LHEAD, user called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
C--------------    End CDE                             --------------
      DIMENSION    LHEADP(9)
      CHARACTER    *(*) CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZDR, 4HOP   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZDROP /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZDROP  ')
#endif

#include "zebra/q_jbyt.inc"


      LHEAD = LHEADP(1)
      IF (LHEAD.EQ.0)        RETURN

#include "zebra/qtrace.inc"
#include "zebra/qstore.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      CALL UOPTC (CHOPT,'LV',IQUEST)
      IFLAG = IQUEST(1)
      IF (IQUEST(2).NE.0)  IFLAG=-1

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LHEAD)
      IF (IQFOUL.NE.0)             GO TO 91
#endif
#if !defined(CERNLIB_QDEBUG)
      IQNS = IQ(KQS+LHEAD-2)
#endif
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.LT.2)             GO TO 19
      WRITE (IQLOG,9018) JQSTOR,LHEAD,IQID,CHOPT
 9018 FORMAT (' MZDROP-  Store',I3,' L/ID=',I9,1X,A4,' Opt=',A)
#endif

   19 KHEAD = LQ(KQS+LHEAD+2)

C--                Drop dependents only, not bank itself

   21 IF   (IFLAG)           22, 31, 41
   22 NS = IQNS
      CALL MZFLAG (IXSTOR,LHEAD,IQDROP,'V')
      CALL VZERO (LQ(KQS+LHEAD-NS),NS)
      GO TO 999

C--                Drop bank + dependents, but not successors

   31 CALL MZFLAG (IXSTOR,LHEAD,IQDROP,'.')
      LN = LQ(KQS+LHEAD)
      IF (LN.EQ.0)                 GO TO 88
      IF (LN.EQ.LHEAD)             GO TO 88
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LN)
      IF (IQFOUL.NE.0)             GO TO 92
#endif
      IF (KHEAD.NE.0)  LQ(KQS+KHEAD)=LN
      LQ(KQS+LN+2)  = KHEAD
      GO TO 999

C--                Drop whole linear structure with all dependents

   41 CALL MZFLAG (IXSTOR,LHEAD,IQDROP,'L')

   88 IF (KHEAD.NE.0)  LQ(KQS+KHEAD)=0

#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

#if defined(CERNLIB_QDEBUG)
   92 NQCASE = 1
      NQFATA = 1
      IQUEST(12) = LN
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(11) = LHEAD
#include "zebra/qtofatal.inc"
#endif
      END
*      ==================================================
#include "zebra/qcardl.inc"
