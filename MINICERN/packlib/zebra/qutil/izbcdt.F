*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:21  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1996/09/23 12:20:34  cernlib
* add Linux
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_Z32)
#include "qutyz32/izbcdt.F"
#elif (defined(CERNLIB_VAXVMS))||(defined(CERNLIB_QMVMI)) \
    ||(defined(CERNLIB_QMDOS))||(defined(CERNLIB_TMO)) \
    ||((defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC)))
#include "qutyinv/izbcdt.F"
#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
#include "not_cern/qutycvx/izbcdt.F"
#elif defined(CERNLIB_CRAY)
#include "not_cern/qutycry/izbcdt.F"
#elif defined(CERNLIB_CDC)
#include "not_cern/qutycdc/izbcdt.F"
#elif defined(CERNLIB_PDP)
#include "not_cern/qutypdp/izbcdt.F"
#else
      SUBROUTINE IZBCDT (NP,ITABT)

C-    This is a very slow stop-gap solution only  !!!

#include "zebra/quest.inc"
#include "zebra/zkrakc.inc"
#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      DIMENSION    NP(9), ITABT(99)


      N     = NP(1)
      LIM   = ITABT(1)
      JGOOD = 0
      JBAD  = 0

      DO 29  JWH=1,N
      JV = IUCOMP (IQHOLK(JWH),IQLETT(1),LIM)

      IF (JV.EQ.0)                 GO TO 27
      JV = ITABT(JV+1)
      IF   (JV+1)            29, 27, 24
   24 JGOOD = JGOOD + 1
      IQCETK(JGOOD) = JV
      GO TO 29

   27 JBAD = JBAD + 1
   29 CONTINUE
      IQUEST(1) = JGOOD
      IQUEST(2) = JBAD
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
