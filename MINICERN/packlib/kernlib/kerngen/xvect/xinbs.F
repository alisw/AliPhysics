*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:34  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMMPW)
#include "mpwgs/xinbs.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/xinbs.F"
#elif defined(CERNLIB_QMIRT)
#include "irtngs/xinbs.F"
#elif defined(CERNLIB_QMIRTD)
#include "irtdgs/xinbs.F"
#else
      SUBROUTINE XINBS (LUN,XAV,NA,XBV,NXB)
C
C CERN PROGLIB# Z203    XINBS           .VERSION KERNFOR  4.31  911111
C ORIG. 15/07/78 JZ
C
      COMMON /SLATE/ NXR, DUMMY(39)
      DIMENSION    XBV(*), XAV(NA)
C
      READ (LUN,IOSTAT=JSTAT) NR,XAV,(XBV(J),J=1,MIN(NR,NXB))
      IF    (JSTAT)          61, 24, 64
   24 NXR = NR
      NR  = MIN(NR,NXB)
      NXB = MAX(1,NR)
      RETURN
C
   61 NXB = 0
      RETURN
C
   64 NXB = -JSTAT
      RETURN
      END
#endif
