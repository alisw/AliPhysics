*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION IBITS (IZW,IZP,NZB)
*
      IBITS=ISHFTR(LSHIFT(IZW,32-IZP-NZB),32-NZB)
*
      RETURN
      END
