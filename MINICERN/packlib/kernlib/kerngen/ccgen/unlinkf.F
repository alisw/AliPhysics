*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      INTEGER FUNCTION UNLINKF (NAME)
C
C CERN PROGLIB# Z265    UNLINKF         .VERSION KERNFOR  4.31  911111
C ORIG. 15/10/91, JZ
C Fortran interface routine to unlink
C
      CHARACTER    NAME*(*)
      INTEGER      UNLINI

      LGNAME  = LNBLNK (NAME)
      UNLINKF = UNLINI (NAME,LGNAME)
      RETURN
      END
