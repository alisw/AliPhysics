*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UFILL (IV,JL,JR,ISTUFF)
C
C CERN PROGLIB# V300    UFILL           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71 JZ
C
      DIMENSION  IV(*)
C
      ITEXT=ISTUFF
         DO 29 J=JL,JR
   29 IV(J)=ITEXT
      RETURN
      END
