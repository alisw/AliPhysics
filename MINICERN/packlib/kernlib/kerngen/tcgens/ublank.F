*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:14  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UBLANK (IV,JL,JR)
C
C CERN PROGLIB# V300    UBLANK          .VERSION KERNFOR  4.08  840613
C ORIG. 01/01/65 JZ
C
      DIMENSION IV(80)
      DATA IBL     /4H     /
C
C
      DO 29 J=JL,JR
   29 IV(J)= IBL
      RETURN
      END
