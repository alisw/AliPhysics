*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVDO1

C-    Subsidiary to ZVERIF, check process table

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcwk.inc"
#include "zebra/zvfwkc.inc"
C--------------    End CDE                             --------------
      PARAMETER   (NVECPR=61)
      DIMENSION    MVECPR(NVECPR)
      EQUIVALENCE (MVECPR,NQSTOR)
      CHARACTER    TEXT*8
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZVDO, 4H1    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZVDO1  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZVDO1   ')
#endif

      IQUEST(1) = LOCF (IQTABV(1)) - 1
      IQUEST(2) = LOCF (LQ(1)) - 1
      IQUEST(3) = LOCF(IQWKTB(1)) - LQASTO
      IQUEST(4) = LOCF(IQWKFZ(1)) - LQASTO
      IQUEST(5) = LOCF(IQDN2(20)) - LQATAB

      NQCASE = 0
      IF (NQSTOR.EQ.-1)            GO TO 71
      IF (NQSTOR.LT.0)             GO TO 81
      IF (IQUEST(1).NE.LQATAB)     GO TO 82
      IF (IQUEST(2).NE.LQASTO)     GO TO 83
      IF (LQBTIS.NE.LQATAB-LQASTO) GO TO 84
      IF (IQUEST(3).NE.LQWKTB)     GO TO 85
      IF (IQUEST(4).NE.LQWKFZ)     GO TO 86
      IF (IQUEST(5).NE.NQTSYS)     GO TO 87
      RETURN

   71 WRITE (IQLOG,9071)
      GO TO 98

   87 NQCASE = 1
   86 NQCASE = NQCASE + 1
   85 NQCASE = NQCASE + 1
   84 NQCASE = NQCASE + 1
   83 NQCASE = NQCASE + 1
   82 NQCASE = NQCASE + 1
   81 NQCASE = NQCASE + 1
      NQFATA = 5

      WRITE (IQLOG,9091)
      TEXT = 'NQSTOR'
      WRITE (IQLOG,9099) TEXT,NQSTOR
      TEXT = 'NQOFFT'
      WRITE (IQLOG,9099) TEXT,NQOFFT
      TEXT = 'NQOFFS'
      WRITE (IQLOG,9099) TEXT,NQOFFS
      TEXT = 'NQALLO'
      WRITE (IQLOG,9099) TEXT,NQALLO
      TEXT = '50->end'
      WRITE (IQLOG,9099) TEXT,(MVECPR(J),J=50,NVECPR)

   98 CONTINUE
#include "zebra/qtraceq.inc"
#include "zebra/qtofatal.inc"

 9091 FORMAT (1X/' ZVDO1.   Process parameters overwritten'
     F/10X,'dump common /MZCA/'/1X)
 9071 FORMAT (1X/' ZVDO1.   No Zebra store initialized.')
#if defined(CERNLIB_B32)
 9099 FORMAT (3X,A,1X,5Z10/(12X,5Z10))
#endif
#if !defined(CERNLIB_B32)
 9099 FORMAT (3X,A,1X,4Z18/(12X,4Z18))
#endif
      END
*      ==================================================
#include "zebra/qcardl.inc"
