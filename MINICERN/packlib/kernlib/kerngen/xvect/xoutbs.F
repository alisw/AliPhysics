*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:35  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIRTD)
#include "irtdgs/xoutbs.F"
#else
      SUBROUTINE XOUTBS (LUN,AV,NA,BV,NB)
C
C CERN PROGLIB# Z203    XOUTBS          .VERSION KERNFOR  3.05  780715
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB), AV(NA)
C
      WRITE (LUN)  NB,AV,BV
      RETURN
      END
#endif
