*
* $Id$
*
* $Log$
* Revision 1.4  1998/09/25 09:33:23  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.3  1997/03/14 17:20:41  mclareni
* WNT mods
*
* Revision 1.2  1996/09/23 12:20:36  cernlib
* add Linux
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:49  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_Z32)
#include "qutyz32/izhnum.F"
#elif (defined(CERNLIB_VAXVMS))||(defined(CERNLIB_DECS)) \
    ||(defined(CERNLIB_DECOSF))||(defined(CERNLIB_MSDOS)) \
    ||(defined(CERNLIB_TMO))   ||((defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))) \
    ||(defined(CERNLIB_WINNT))
#include "qutyinv/izhnum.F"
#elif (defined(CERNLIB_CONVEX))||(defined(CERNLIB_SINGLE))
#include "not_cern/qutycvx/izhnum.F"
#elif defined(CERNLIB_CRAY)
#include "not_cern/qutycry/izhnum.F"
#elif defined(CERNLIB_CDC)
#include "not_cern/qutycdc/izhnum.F"
#elif defined(CERNLIB_PDP)
#include "not_cern/qutypdp/izhnum.F"
#endif
