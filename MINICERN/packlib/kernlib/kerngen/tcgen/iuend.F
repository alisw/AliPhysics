*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUEND (NDA)
C
C CERN PROGLIB# M410    IUEND           .VERSION KERNFOR  1.0   670527
C ORIG. 27/05/67
C
C-    RETURNS THE SEARCH-PARAMETERS IF IUFORW,IUBACK,ULEFT,URIGHT
C
      COMMON /SLATE/ ND,JEND, DUMMY(38)
C
      NDA  = ND
      IUEND= JEND
      RETURN
      END
