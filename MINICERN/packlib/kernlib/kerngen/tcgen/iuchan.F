*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUCHAN (X,XL,DX,NX)
C
C CERN PROGLIB# Y201    IUCHAN          .VERSION KERNFOR  3.15  820113
C ORIG. 01/01/66 JZ
C
#if defined(CERNLIB_B36M)
      AN= (X-XL)/DX + 1.000001
#endif
#if !defined(CERNLIB_B36M)
      AN= (X-XL)/DX + 1.00001
#endif
      AN= MAX(AN,0.)
      AN= MIN(AN, FLOAT(NX+1))
      N=AN
      IF (N.LT.0) N=NX+1
      IUCHAN=N
      RETURN
      END
