*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/chdirf.F"
#else
      INTEGER FUNCTION CHDIRF (NAME)
C
C CERN PROGLIB# Z265    CHDIRF          .VERSION KERNFOR  4.26  910313
C ORIG. 22/02/91, JZ
C Fortran interface routine to chdir
C
      CHARACTER    NAME*(*)
      INTEGER      CHDIRI

      LGNAME = LNBLNK (NAME)
      CHDIRF = CHDIRI (NAME,LGNAME)
      RETURN
      END
#endif
