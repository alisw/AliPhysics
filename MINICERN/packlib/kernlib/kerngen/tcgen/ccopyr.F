*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CCOPYR (FR,TO,NCH)
C
C CERN PROGLIB# M432    CCOPYR          .VERSION KERNFOR  4.39  940228
C ORIG. 08/10/91, JZ
C
C-    Copy NCH characters from FR to TO shifting right

      CHARACTER    FR(512)*1, TO(512)*1

      DO 24  J=NCH,1,-1
   24 TO(J) = FR(J)
      RETURN
      END
