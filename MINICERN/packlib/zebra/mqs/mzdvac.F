*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:00  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      FUNCTION MZDVAC (IXDIVP)

C-    collect and return bits of divisions requested by IXDIVP,
C-    system called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXDIVP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZDV, 4HAC   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZDVAC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZDVAC  ')
#endif

#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit1.inc"
#include "zebra/q_jbytet.inc"


      IXIN = IXDIVP(1)

#include "zebra/qtrace.inc"

      JST = JBYT (IXIN,27,6)
      IF    (JST.EQ.JQSTOR)        GO TO 31
      IF (JST-16.EQ.JQSTOR)        GO TO 21
      CALL MZSDIV (IXIN,-7)
      IF (JST.LT.16)               GO TO 31

C--                COMPOSITE INDEX

   21 IXAC = JBYT (IXIN,1,20)
      IXGE = JBYT (IXIN,21,6)
      IF (IXGE.EQ.0)               GO TO 59
      IF (IXGE.LT.16)              GO TO 41

C--                ERROR TERMINATION

   29 CALL MZSDIV (IXIN,0)

C--                SINGLE DIVISION INDEX

   31 JDIV = JBYT (IXIN,1,26)
      IF (JDIV.GE.25)              GO TO 29
      IXAC = 0
      IF (JDIV.GE.21)              GO TO 33
      IXAC = MSBIT1 (IXAC,JDIV)
      GO TO 59

   33 IXGE = MSBIT1 (0, JDIV-20)

C----              COLLECT BITS FOR GENERIC INDEX

   41 JDIV = 1
   42 IF (JDIV.EQ.JQDVLL+1)  JDIV=JQDVSY
      IF (JBYTET(IXGE,IQKIND(KQT+JDIV),21,4).EQ.0)  GO TO 47
      IXAC = MSBIT1 (IXAC,JDIV)
   47 JDIV = JDIV + 1
      IF (JDIV.LT.21)              GO TO 42

   59 MZDVAC = IXAC
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.LT.11)            GO TO 999
      WRITE (IQLOG,9859) IXIN, IXAC
#endif
#if (defined(CERNLIB_QDEVZE))&&(defined(CERNLIB_HEX))
 9859 FORMAT (1X/' DEVZE MZDVAC,  IXDIVP,MZDVAC=',2(2X,Z8))
#endif
#if (defined(CERNLIB_QDEVZE))&&(!defined(CERNLIB_HEX))
 9859 FORMAT (1X/' DEVZE MZDVAC,  IXDIVP,MZDVAC=',2(2X,O11))
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
