*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CROSS(X,Y,Z)
C
C CERN PROGLIB# F117    CROSS           .VERSION KERNFOR  2.02  730125
C ORIG. 01/01/64
C
      COMMON /SLATE/Z1,Z2,XQX(38)
      DIMENSION X(3),Y(3),Z(3)
C
C
      Z1=X(2)*Y(3)-X(3)*Y(2)
      Z2=X(3)*Y(1)-X(1)*Y(3)
      Z(3)=X(1)*Y(2)-X(2)*Y(1)
      Z(1)=Z1
      Z(2)=Z2
      RETURN
      END
