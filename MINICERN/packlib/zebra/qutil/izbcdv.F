*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:22  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1996/09/23 12:20:35  cernlib
* add Linux
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
*       .
#include "zebra/pilot.h"
#if defined(CERNLIB_Z32)
#include "qutyz32/izbcdv.F"
#elif (defined(CERNLIB_VAXVMS))||(defined(CERNLIB_QMVMI)) \
    ||(defined(CERNLIB_QMDOS))||(defined(CERNLIB_TMO)) \
	||((defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC)))
#include "qutyinv/izbcdv.F"
#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_SINGLE))
#include "not_cern/qutycvx/izbcdv.F"
#elif defined(CERNLIB_CRAY)
#include "not_cern/qutycry/izbcdv.F"
#elif defined(CERNLIB_CDC)
#include "not_cern/qutycdc/izbcdv.F"
#elif defined(CERNLIB_PDP)
#include "not_cern/qutypdp/izbcdv.F"
#else
      SUBROUTINE IZBCDV (NP)

C-    This is a very slow stop-gap solution only  !!!

#include "zebra/zkrakc.inc"
#include "zebra/zbcd.inc"
C--------------    END CDE                             --------------
      DIMENSION    NP(9)


      N = NP(1)

      DO 29  JWH=1,N
      JV = IUCOMP (IQHOLK(JWH),IQLETT(1),63)
      IF (JV.EQ.0)  JV=57
   29 IQCETK(JWH) = JV
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
