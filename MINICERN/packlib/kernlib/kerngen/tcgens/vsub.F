*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VSUB (A,B,X,N)
C
C CERN PROGLIB# F121    VSUB            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),B(*),X(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 X(I)= A(I) - B(I)
      RETURN
      END
