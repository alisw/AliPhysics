*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION AMAXMU (A,IDO,IW,NA)
C
C CERN PROGLIB# E103    AMAXMU          .VERSION KERNFOR  1.0   680315
C ORIG. 15/03/68
C
      DIMENSION A(99)
C
      AM=0.
      II=1
C
         DO 19 JN=1,IDO
      I=II
C
         DO 17 JK=1,IW
      AM = MAX (AM, ABS(A(I)))
   17 I=I+1
   19 II=II+NA
C
      AMAXMU=AM
      RETURN
      END
