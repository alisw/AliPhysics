*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CFWEOF (LUNPTR, MEDIUM, NEOF)
C CERN PROGLIB#         CFWEOF          .VERSION KERNFOR  4.26  910313
C     write NEOF system file marks, dummy routine
C      LUNPTR  file pointer
C      MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
C     *NEOF*   number of file marke to be writtin / actually written

      NEOF = 0
      RETURN
      END
