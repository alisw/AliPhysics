*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VBLANK (A,N)
C
C CERN PROGLIB# F121    VBLANK          .VERSION KERNFOR  4.08  840613
C ORIG. 01/07/71
C
      INTEGER      A(*), BLANK
      DATA  BLANK  / 4H     /
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= BLANK
      RETURN
      END
