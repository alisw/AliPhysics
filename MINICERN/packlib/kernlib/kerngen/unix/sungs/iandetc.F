*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:08  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IAND (I,J)
C
C CERN PROGLIB# M441    IANDETC         .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
C     SUN versions of functions  IAND, IOR, IEOR

      IAND = AND(I,J)
      END
      FUNCTION IOR (I,J)

      IOR = OR(I,J)
      END
      FUNCTION IEOR (I,J)

      IEOR = XOR(I,J)
      END
