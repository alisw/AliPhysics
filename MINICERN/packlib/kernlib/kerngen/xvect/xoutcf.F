*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:35  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE XOUTCF (LUN,BV,NB)
C
C CERN PROGLIB# Z203    XOUTCF          .VERSION KERNFOR  3.15  820113
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB)
C
      WRITE (LUN,8000) BV
      RETURN
C
#if defined(CERNLIB_A4)
 8000 FORMAT (35A4)
#endif
#if defined(CERNLIB_A5)
 8000 FORMAT (28A5)
#endif
#if defined(CERNLIB_A6)
 8000 FORMAT (23A6)
#endif
#if defined(CERNLIB_A8)
 8000 FORMAT (16A8)
#endif
#if defined(CERNLIB_A10)
 8000 FORMAT (14A10)
#endif
      END
