*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_TCGEN)
      SUBROUTINE POWEZE (X,A,N)
C
C CERN PROGLIB# B400    POWEZE          .VERSION KERNFOR  2.09  751101
C ORIG. 20/08/75
C
      DIMENSION A(2)
      A(1)= 1.
      IF (N.LE.1) RETURN
      DO 10 I= 2,N
   10 A(I)= X * A(I-1)
      RETURN
      END
#endif
