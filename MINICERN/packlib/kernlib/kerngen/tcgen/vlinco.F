*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:56  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VLINCO (B,ALPHA,C,BETA,A,N)
C
C CERN PROGLIB# F121    VLINCO          .VERSION KERNFOR  4.12  850702
C ORIG. 01/07/71 VLINE; CHANGE NAME TO VLINCO  02/07/85
C
      DIMENSION A(*),B(*),C(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= B(I)*ALPHA + C(I)*BETA
C
      RETURN
      END
