*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VFIX (A,IA,N)
C
C CERN PROGLIB# F121    VFIX            .VERSION KERNFOR  3.15  820113
C ORIG. 03/05/72
C
      DIMENSION A(*), IA(*)
C
      IF (N.LE.0)  RETURN
      DO 9 J=1,N
    9 IA(J)= A(J)
      RETURN
      END
