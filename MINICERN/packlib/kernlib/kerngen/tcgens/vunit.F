*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VUNIT (A,B,N)
C
C CERN PROGLIB# F121    VUNIT           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),B(*)
C
      IF (N.LE.0)  RETURN
      XX= 0.
      DO 9 I=1,N
    9 XX= XX + A(I)*A(I)
      IF (XX.LE.0.)  RETURN
      YY=1./ SQRT (XX)
      DO 19 I= 1,N
   19 B(I)=A(I) * YY
      RETURN
      END
