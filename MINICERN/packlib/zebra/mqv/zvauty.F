*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVAUTY

C-    Execute automatic verification

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/zvfaut.inc"
#include "zebra/zvfwkc.inc"
C--------------    END CDE                             --------------
#if !defined(CERNLIB_QTRHOLL)||defined(CERNLIB_A6M)
      PARAMETER    (LEVSEC = 2)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      PARAMETER    (LEVSEC = 4)
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZVAU, 4HTY   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZVAUTY /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZVAUTY  ')
#endif

      IF (IQVSTA.EQ.0)             RETURN

C--                Check this system call is a primary call

      IF (NQTRAC.GE.LEVSEC)        RETURN

C--                Check primary store exists

      IF (NQSTOR.EQ.-1)            RETURN

C--                Step the Verify ID

      IF (IQVID(2).EQ.2147483647)  THEN
          IQVID(1) = IQVID(1) + 1
          IQVID(2) = 0
        ENDIF

      IQVID(2) = IQVID(2) + 1

      IF (IQVSTA.GE.0)             GO TO 41

C--                Running in dummy mode

      IF (IQVID(1).NE.IQVTHR(1))   RETURN
      IF (IQVID(2).NE.IQVTHR(2))   RETURN

C--                Switch to active mode

      IQVSTA = 1
      IF (IQVLOG.EQ.0)             GO TO 41

      DO 24  JST=1,NQSTOR+1
      K = NQOFFT(JST)
   24 IQTABV(K+10) = IQVLOG

C--                Running in active mode

   41 JSTO = JQSTOR
      JDIV = JQDIVI

#include "zebra/qtraceq.inc"
      CALL ZVERIF (-2, 0, 'Automatic verification')
#include "zebra/qtrace99.inc"

      IF (JQSTOR.NE.JSTO)  THEN
          IXSTOR = 0
          CALL SBYT (JSTO,IXSTOR,27,4)
          CALL MZSDIV (IXSTOR,-7)
        ENDIF
      JQDIVI = JDIV
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
