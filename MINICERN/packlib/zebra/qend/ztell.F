*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZTELL (IDP,JFLP)

C-    Signal trouble situation,
C-    called from the garbage collector if not enough space
C-    may be user called with IDP > 100

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/quest.inc"
#include "zebra/mzca.inc"
#include "zebra/mzcb.inc"
#include "zebra/fzcx.inc"
C--------------    END CDE                             --------------
      COMMON /ZTELLC/  ID,JFLUS

      DIMENSION    IDP(9), JFLP(9)

      PARAMETER   (NLIST = 3 )
      DIMENSION    LIST(NLIST), JFLDEF(NLIST)
#if defined(CERNLIB_QMVDS)
      SAVE         LIST, JFLDEF
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZTEL, 4HL    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZTELL  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZTELL   ')
#endif
      DATA  LIST   / 19, 61, 62 /
      DATA  JFLDEF /  2,  2,  2 /

C-    Return code:   0 return, 1 QNEXT, 2 ZFATAL, 3 ZEND


      ID    = IDP(1)
      JFLUS = JFLP(1)
      JFLSV = MAX (JFLUS,0)
      IQUEST(1) = ID
#if defined(CERNLIB_QPRINT)
      LOGLEV = NQLOGD
      IF (ID.GE.11)  THEN
          IF (ID.LE.19)  LOGLEV = MAX (LOGLEV,LOGLVX)
          IF (ID.EQ.15)  LOGLEV = MAX (LOGLEV,NQLOGL)
          IF (ID.EQ.99)  LOGLEV = MAX (LOGLEV,NQLOGL)
        ENDIF
      IF (LOGLEV.GE.1)   WRITE (IQLOG,9001)  ID,JFLUS
 9001 FORMAT (' ZTELL.  ----   Called for Condition',I5,',',I4)
#endif

C--                System conditions

      IF (ID.GE.100)               GO TO 31
      J = IUCOMP (ID,LIST,NLIST)
      IF (J.NE.0)   THEN
          JFLUS = JFLDEF(J)
        ENDIF

C--                CALL ZTELUS

   31 CALL ZTELUS
      JFLUS = MAX (JFLUS,JFLSV)
      JFLUS = MIN (JFLUS,3) + 1
      GO TO (40,41,42,43), JFLUS

   41 IF (NQPHAS.LE.0)             GO TO 42
      NQTRAC = 0
      CALL QNEXTE

   43 NQPHAS = -2
      CALL ZEND

   42 CONTINUE
#include "zebra/qtrace.inc"
      CALL ZFATAM ('EXIT VIA ZTELL.')

   40 RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
