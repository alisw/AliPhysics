*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:16  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE JZWIND (CHPA1)

C-    Reset processor for trap transfer

#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/jzuc.inc"
#include "zebra/jzc.inc"
C--------------    END CDE                             -----------------  ------
      CHARACTER    CHPA1*4
#if defined(CERNLIB_A4)
      CHARACTER    CHNMTO*4, CHIAM*4
#endif
#if defined(CERNLIB_A8)
      CHARACTER    CHNMTO*8, CHIAM*8
#endif
#if defined(CERNLIB_EQUHOLCH)
      INTEGER      IAMTO
      EQUIVALENCE (CHNMTO, IAMTO), (CHIAM, NQME(1))
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HJZWI, 4HND   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HJZWIND /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'JZWIND  ')
#endif

#if defined(CERNLIB_QDEBUG)
#include "zebra/q_jbit.inc"
#endif

      IF (JQLEV.LE.0)              RETURN

#include "zebra/qtraceq.inc"

      CHNMTO  = CHPA1
#if defined(CERNLIB_EQUHOLCH)
   21 IF (NQME(1).NE.IAMTO)  THEN
          IF (JQLEV.EQ.0)          GO TO 91
          CALL JZOUT ('??? ')
          GO TO 21
        ENDIF
#endif
#if !defined(CERNLIB_EQUHOLCH)
   21 CALL UHTOC (NQME(1),4,CHIAM,4)
      IF (CHIAM(1:4).NE.CHNMTO(1:4))  THEN
          IF (JQLEV.EQ.0)          GO TO 91
          CALL JZOUT ('??? ')
          GO TO 21
        ENDIF
#endif

#if defined(CERNLIB_QDEBUG)
      IF (JBIT(JQLLEV,13).NE.0)  WRITE (IQLOG,9016) CHNMTO(1:4)
 9016 FORMAT    (/' =======  JZWIND to "',A,' done.')
#endif
      IQUEST(1) = 0
#include "zebra/qtrace99.inc"
      RETURN

   91 IQUEST(1) = 7
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
