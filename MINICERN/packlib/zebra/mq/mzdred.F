*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:19  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZDRED (IXDIVP)

C-    Reduce space reserved for division to initial size,
C-    but not smaller than current occupation

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIVP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZDR, 4HED   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZDRED /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZDRED  ')
#endif


#include "zebra/qtrace.inc"

      CALL MZSDIV (IXDIVP,4)
      IF (JQDIVI.LT.3)             GO TO 999
      CALL MZRESV

#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.1)  WRITE (IQLOG,9028) JQSTOR,JQDIVI,NQRESV
 9028 FORMAT (' MZDRED-  called for Store/Div',2I3,' Free',I7)
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)  WRITE (IQLOG,9829) JQDIVI,JQSHAR,NQRESV
 9829 FORMAT (1X/' DEVZE MZDRED entered,  JQDIVI,JQSHAR,NQRESV= ',3I8)
#endif
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif
C----              Construct Memory Occupation table

      MQDVGA = 0
      MQDVWI = 0
      JQSTMV = -1
      CALL MZTABM

C----              Shift divisions

      IF (JQSHAR.NE.0)  THEN
          NFREE = MIN (NQRESV,
     +                 LQEND(KQT+JQSHR2)-LQSTA(KQT+JQSHR1)
     +                 - (NQDINI(KQT+JQDIVI)+NQDINI(KQT+JQDIVN)) )
        ELSE
          NFREE = NQRESV - MAX (0, NQDINI(KQT+JQDIVI) -
     +                   (LQEND(KQT+JQDIVI)-LQSTA(KQT+JQDIVI)) )
        ENDIF

      IF (NFREE.LE.24)             GO TO 999

C--                Forward division

      IF (JQMODE.EQ.0)  THEN
          JQDVM2 = JQDIVI
        ELSE

C--                Reverse division

          JQDVM2 = JQDIVN
        ENDIF

      JQDVM1 = 2
      JQSTMV = JQSTOR
      NQDVMV = NFREE

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)  WRITE (IQLOG,9836) JQDVM1,JQDVM2,NQDVMV
 9836 FORMAT (' DEVZE MZDRED,  JQDVM1,JQDVM2,NQDVMV=    ',3I8)
#endif
C--                Construct Link Relocation table

      CALL MZTABR
      CALL MZTABS

C----              Relocate + memory move

      CALL MZTABX
      CALL MZTABF
      CALL MZRELX
      CALL MZMOVE
      NQDRED(KQT+JQDIVI) = NQDRED(KQT+JQDIVI) + 1
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
