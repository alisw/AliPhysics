*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:41  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE BLOW (CHAR,INT,NCHP,NBITS)
C
C CERN PROGLIB# M426    BLOW            .VERSION KERNFOR  4.31  911111
C
      DIMENSION    CHAR(*), INT(*)
      DIMENSION    NCHP(9), IPAR(4)
      DATA  IPAR   / 0, 0, 0, 0 /
C
      NCH     = NCHP(1)
      IPAR(1) = NBITS
      IPAR(3) = NCH*IPAR(1) + 127
C
      CALL UPKCH (CHAR,INT,NCH,IPAR)
      RETURN
      END
