*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION PROXIM (ANGL,STAND)
C
C CERN PROGLIB# V306    PROXIM          .VERSION KERNFOR  4.21  890323
C ORIG. 15/03/68 JZ, re-write 8/03/89 K.S.Koelbig
C
      PARAMETER (PI2 = 6.28318 53071 79586D0, RPI2 =1/PI2)

      PROXIM = ANGL + PI2*ANINT(RPI2*(STAND-ANGL))
      RETURN
      END
