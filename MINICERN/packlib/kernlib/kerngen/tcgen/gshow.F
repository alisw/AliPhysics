*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:41  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxgs/gshow.F"
#else
      SUBROUTINE GSHOW (LUN)
C
C CERN PROGLIB# Z001    GSHOW           .VERSION KERNFOR  4.29  910718
C ORIG. 18/12/78
C
      CALL KERNGT (LUN)
      RETURN
      END
#endif
