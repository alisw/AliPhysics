*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE  VCOPYN (A,X,N)
C
C CERN PROGLIB# F121    VCOPYN          .VERSION KERNFOR  3.15  820113
C ORIG. 01/04/73
C
      DIMENSION A(*), X(*)
C
      IF(N.LE.0)  RETURN
      DO 10  I=1,N
   10 X(I)= -A(I)
      RETURN
      END
