*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:25  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1996/09/23 12:20:38  cernlib
* add Linux
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_Z32)
#include "qutyz32/zitoh.F"
#elif (defined(CERNLIB_VAXVMS))||(defined(CERNLIB_QMVMI)) \
    ||(defined(CERNLIB_QMDOS))||(defined(CERNLIB_TMO)) \
    ||((defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC)))
#include "qutyinv/zitoh.F"
#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
#include "not_cern/qutycvx/zitoh.F"
#elif defined(CERNLIB_CRAY)
#include "not_cern/qutycry/zitoh.F"
#elif defined(CERNLIB_CDC)
#include "not_cern/qutycdc/zitoh.F"
#elif defined(CERNLIB_PDP)
#include "not_cern/qutypdp/zitoh.F"
#else
      SUBROUTINE ZITOH (INTV,HOLL,NP)

C-    This is a rather slow stop-gap solution only  !!!

#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      COMMON /SLATE/ DUMMY(8), MM(4), DUMB(28)

      DIMENSION    INTV(99), HOLL(99), NP(9)
      DIMENSION    MPAK(2)
      DATA  MPAK   /6,4/


      N = NP(1)

      DO 39  JW=1,N
      CALL UPKBYT (INTV(JW),1,MM(1),4,MPAK(1))

      DO 16  J=1,4
      JV = MM(J)
      IF (JV.EQ.0)  JV=45
   16 MM(J) = IQLETT(JV)

      CALL UBUNCH (MM(1),HOLL(JW),4)
   39 CONTINUE
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
