*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:41  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE GSHOWF (LUNP)
C
C CERN PROGLIB# Z001    GSHOWF          .VERSION KERNFOR  4.39  940228
C ORIG. 18/12/78
C
      COMMON /SLATE/ ICH(32),DUM(8)
      CHARACTER    VIDQQ*(*)
      PARAMETER   (VIDQQ = '@(#)' //
     + ' KERNFOR  4.42  1995/10/11 17.43        >')
      PARAMETER    (LUNSTP=6)

      LUN = LUNP
      IF (LUN.EQ.0)  LUN = LUNSTP

      WRITE (LUN,9002) ICH
      CALL CLEFT (VIDQQ,1,0)
      RETURN

 9002 FORMAT (1X/25H KERNGT.  KERNLIB from:  , 32A1/25X,
     + ' KERNFOR  4.42  1995/10/11 17.43        ')
      END
