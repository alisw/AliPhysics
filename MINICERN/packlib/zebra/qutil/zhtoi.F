*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:24  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1996/09/23 12:20:37  cernlib
* add Linux
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_Z32)
#include "qutyz32/zhtoi.F"
#elif (defined(CERNLIB_VAXVMS))||(defined(CERNLIB_QMVMI)) \
    ||(defined(CERNLIB_QMDOS))||(defined(CERNLIB_TMO)) \
    ||((defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC)))
#include "qutyinv/zhtoi.F"
#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
#include "not_cern/qutycvx/zhtoi.F"
#elif defined(CERNLIB_CRAY)
#include "not_cern/qutycry/zhtoi.F"
#elif defined(CERNLIB_CDC)
#include "not_cern/qutycdc/zhtoi.F"
#elif defined(CERNLIB_PDP)
#include "not_cern/qutypdp/zhtoi.F"
#else
      SUBROUTINE ZHTOI (HOLL,INTV,NP)

C-    This is a slow stop-gap solution only  !!!

#include "zebra/zceta.inc"
C--------------    END CDE                             --------------
      COMMON /SLATE/ DUMMY(8), MM(4), MIM(4), DUMB(24)

      DIMENSION    HOLL(99), INTV(99), NP(9)
      DIMENSION    MPAK(2)
      DATA  MPAK   /6,4/


      N = NP(1)

      DO 39  JW=1,N
      CALL UBLOW  (HOLL(JW),MM(1),4)
      CALL IZHNUM (MM,MIM,4)
      INTW = 0

      MM(1) = IQTCET(MIM(1))

      DO 16  J=2,4
      MM(J) = IQTCET(MIM(J))
      IF (MM(J).EQ.45)  MM(J) = 0
   16 CONTINUE

      CALL PKBYT (MM(1),INTW,1,4,MPAK(1))
   39 INTV(JW) = INTW
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
