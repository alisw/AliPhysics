*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:55  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZNEED (IXDIV,NEEDP,CHOPT)

C-    Check that NEEDP words are available in the specified division,
C-    without garbage collection, unless option G is given;
C-    return in IQUEST(11) the space available beyond this,
C-    i.e. IQUEST(11) negative if not enough space;
C-    return IQUEST(12/13) = current/maximum size of division;
C-    user called

#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXDIV(9),NEEDP(9)
      CHARACTER    *(*) CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZNE, 4HED   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZNEED /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZNEED  ')
#endif

#include "zebra/q_jbyt.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      JDV  = IXDIV(1)
      NEED = NEEDP(1)
      IF (JBYT(JDV,27,6).NE.JQSTOR)  GO TO 22
      JQDIVI = JBYT (JDV,1,26)
      IF (JQDIVI.EQ.0)               GO TO 22
      IF (JQDIVI.LT.21)              GO TO 23
   22 CALL MZSDIV (JDV,4)
   23 CONTINUE

C--            Allocate space

      CALL MZRESV
      NQRESV = NQRESV - NEED
      IF (NQRESV.LT.0)             GO TO 41

   28 IQUEST(11) = NQRESV
      IQUEST(12) = LQEND(KQT+JQDIVI) - LQSTA(KQT+JQDIVI)
      IQUEST(13) = NQDMAX(KQT+JQDIVI)
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.2)
     + WRITE (IQLOG,9029) JQSTOR,JQDIVI,NEED,NQRESV,CHOPT
 9029 FORMAT (' MZNEED-  Store/Div',2I3,' NEED/Excess=',2I8
     F,' Opt=',A)
#endif
#include "zebra/qtrace99.inc"
      RETURN

C----              Garbage collection

   41 CALL UOPTC (CHOPT,'G',IQUEST)
      IF (IQUEST(1).EQ.0)          GO TO 28
      NQPERM = 1
      CALL MZGAR1
      NQPERM = 0
      GO TO 28
      END
*      ==================================================
#include "zebra/qcardl.inc"
