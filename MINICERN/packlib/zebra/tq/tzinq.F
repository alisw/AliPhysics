*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:14:54  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZINQ  (IXSTOR,IXTITL,LTT,IFLAG)

C-    Inquire about the title structure in store IXSTOR
C-        return:  the index of the title division
C-                 the adr of the last (IFLAG=0) or the first bank

#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
C--------------    END CDE                             -----------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HTZIN, 4HQ    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HTZINQ  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'TZINQ   ')
#endif

#include "zebra/q_jbyt.inc"
#include "zebra/q_sbyt.inc"

#include "zebra/qtraceq.inc"

#include "zebra/qstore.inc"

      L = LQT(KQT+1)
      IF (L.EQ.0)                  GO TO 70

      JDV  = MZFDIV (-7, L)
      IF (JDV.EQ.JQDVSY)  JDV = 24
      IXDV = MSBYT (JQSTOR,JDV,27,4)
      IF (IFLAG.NE.0)              GO TO 71

C--       go to the end

   24 IF (LQ(KQS+L).EQ.0)             GO TO 71
      L = LQ(KQS+L)
      GO TO 24

   70 IXDV = MSBYT (JQSTOR,0,27,4)
   71 IXTITL = IXDV
      LTT    = L

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
