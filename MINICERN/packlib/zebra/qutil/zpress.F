*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:43  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZPRESS (IXSTOR,LGOP)

C-    REMOVE DEAD BANKS FROM (CIRCULAR) LINEAR STRUCTURE

#include "zebra/mqsysh.inc"
C--------------    END CDE                             --------------
      DIMENSION    LGOP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZPRE, 4HSS   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZPRESS /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZPRESS  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      LGO = LGOP(1)
      IF (LGO.EQ.0)          RETURN

#include "zebra/qtraceq.inc"
#include "zebra/qstore.inc"

      L   = LGO
      KGO = LQ(KQS+L+2)
      K   = KGO

   11 IF (JBIT(IQ(KQS+L),IQDROP).NE.0)  GO TO 21
   12 K = L
      L = LQ(KQS+K)
      IF (L.EQ.0)                  GO TO 999
      IF (L.NE.LGO)                GO TO 11
      LQ(KQS+K) = LQ(KQS+KGO)
#include "zebra/qtrace99.inc"
      RETURN

C--                SERIES OF DEAD BANKS

   21 L = LQ(KQS+L)
      IF (L.EQ.0)                  GO TO 24
      IF (L.EQ.LGO)                GO TO 25
      IF (JBIT(IQ(KQS+L),IQDROP).NE.0)  GO TO 21
      LQ(KQS+K)   = L
      LQ(KQS+L+2) = K
      GO TO 12

   24 LQ(KQS+K) = 0
      GO TO 999

   25 LQ(KQS+K) = LQ(KQS+KGO)
      IF (K.NE.KGO)                GO TO 999
      LQ(KQS+KGO) = 0
      GO TO 999
      END
*      ==================================================
#include "zebra/qcardl.inc"
