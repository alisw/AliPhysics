*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:14:55  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZSHUN (IXSTOR,LBK,IFLAG)

C-    Shunt new banks(s) into the title structure

#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
C--------------    END CDE                             -----------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HTZSH, 4HUN   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HTZSHUN /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'TZSHUN  ')
#endif

#include "zebra/q_jbyt.inc"

#include "zebra/qtraceq.inc"

#include "zebra/qstore.inc"

      L = LQT(KQT+1)
      IF (L.EQ.0)                  GO TO 31
      IF (IFLAG.NE.0)              GO TO 31

C--       go to the end

   24 IF (LQ(KQS+L).NE.0)  THEN
          L = LQ(KQS+L)
          GO TO 24
        ENDIF

      CALL ZSHUNT (IXSTOR,LBK,L,0,1)
      GO TO 999

C--       connect at the start

   31 CALL ZSHUNT (IXSTOR,LBK,LQT(KQT+1),1,1)

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
