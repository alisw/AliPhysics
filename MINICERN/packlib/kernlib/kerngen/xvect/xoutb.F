*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:34  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIRTD)
#include "irtdgs/xoutb.F"
#else
      SUBROUTINE XOUTB  (LUN,BV,NB)
C
C CERN PROGLIB# Z203    XOUTB           .VERSION KERNFOR  3.05  780715
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB)
C
      WRITE (LUN)  NB,BV
      RETURN
      END
#endif
