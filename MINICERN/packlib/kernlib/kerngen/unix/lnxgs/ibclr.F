*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:07  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IBCLR(IZW,IZP)
      IBCLR = AND (IZW,NOT(LSHIFT(1,IZP)))
      RETURN
      END
