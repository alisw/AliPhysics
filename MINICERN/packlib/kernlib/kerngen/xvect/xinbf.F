*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:34  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIRT)
#include "irtngs/xinbf.F"
#elif defined(CERNLIB_QMIRTD)
#include "irtdgs/xinbf.F"
#else
      SUBROUTINE XINBF (LUN,XBV,NXB)
C
C CERN PROGLIB# Z203    XINBF           .VERSION KERNFOR  4.31  911111
C ORIG. 15/07/78 JZ
C
      DIMENSION    XBV(NXB)
C
      READ (LUN,IOSTAT=JSTAT) XBV
      IF    (JSTAT)          61, 24, 64
C
   61 NXB = 0
      RETURN
C
   64 NXB = -JSTAT
   24 RETURN
      END
#endif
