*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:07  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
        FUNCTION IRNDM(ISEEK)
        COMMON /MCGN/ MCGN

C--       Integer random number
      MCGN  = MCGN * 69069
      IRNDM = ISHFT (MCGN,-1)
      RETURN
        END
