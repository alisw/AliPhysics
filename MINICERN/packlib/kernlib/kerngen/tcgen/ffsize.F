*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMHPX)
#include "hpxgs/ffsize.F"
#elif defined(CERNLIB_QMSGI)
#include "sgigs/ffsize.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/ffsize.F"
#else
      SUBROUTINE FFSIZE (LUN,NWREC,JREC,ISTAT)
C
C CERN PROGLIB#         FFSIZE  DUMMY   .VERSION KERNFOR  4.26  910313
C
C-    Get the size of a disk file connected to Fortran unit LUN
C-    dummy version


      ISTAT = -1
      RETURN
      END
#endif
