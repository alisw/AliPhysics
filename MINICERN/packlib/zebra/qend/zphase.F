*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZPHASE (NPHP)

C-    CHANGE PROGRAM PHASE,  USER CALLED

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    NPHP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZPHA, 4HSE   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZPHASE /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZPHASE  ')
#endif


      NPH = NPHP(1)
      IF (NPH.EQ.0)  NPH=1
      IF   (NQPHAS)          41,31,21

C--                CURRENT STATE IS OPERATION

   21 IF (NPH.LT.0)                GO TO 51
      NQPHAS = NPH
      RETURN

C--                CURRENT STATE IS INIT

   31 IF (NPH.LT.0)                GO TO 51
      NQPHAS = NPH
#if defined(CERNLIB_QPRINT)
      IF (NQLOGD.LT.-1)      RETURN
      WRITE (IQLOG,9032) NPH
 9032 FORMAT (1X/' ZPHASE.  Start Operation Phase',I5)
#endif
      RETURN

C--                CURRENT STATE IS TERM

   41 IF (NPH.LT.0)                GO TO 52
#if defined(CERNLIB_QPRINT)
      IF (NQLOGD.LT.-2)      RETURN
      WRITE (IQLOG,9042) NPH
 9042 FORMAT (1X/' !!!!! ZPHASE refuses to go back from Termination'
     F,' to Operation Phase',I5)
#endif
      RETURN

C----              START TERMINATION PHASE

#if !defined(CERNLIB_QPRINT)
   51 CONTINUE
   52 CONTINUE
#endif
#if defined(CERNLIB_QPRINT)
   51 IF (NQLOGD.GE.-1)  WRITE (IQLOG,9051)
 9051 FORMAT (1X/' ZPHASE.  Start Termination Phase.')

   52 IF (NQLOGD.GE.-1)  WRITE (IQLOG,9052) NPH
 9052 FORMAT (1X/' ZPHASE.  Termination mode',I5)
#endif
      NQPHAS = -1

#include "zebra/qtrace.inc"

      CALL MZWORK (0,0,0,-1)
      IF (NPH.GE.-1)               GO TO 999

      CALL MZGARB (JQDVSY,21)
      IF (NPH.GE.-2)               GO TO 999

      LA = LQSTA(JQDVSY)

      JDIV = 1
   56 JDIV = JDIV + 1
      LQSTA(JDIV) = LA
      LQEND(JDIV) = LA
      NQDRED(JDIV) = NQDRED(JDIV) + 1
      IF (JDIV.LT.JQDVLL)          GO TO 56

      IX = MZIXCO (22,23,0,0)
      CALL MZGARB (IX,0)
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
