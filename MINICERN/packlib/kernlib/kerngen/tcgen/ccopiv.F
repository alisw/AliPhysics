*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CCOPIV (FR,TO,NCH)
C
C CERN PROGLIB# M432    CCOPIV          .VERSION KERNFOR  4.39  940228
C ORIG. 20/01/93, JZ
C
C-    Copy NCH characters from FR to TO inverted

      CHARACTER    FR(512)*1, TO(512)*1

      DO 24  J=1,NCH
   24 TO(J) = FR(NCH+1-J)
      RETURN
      END
