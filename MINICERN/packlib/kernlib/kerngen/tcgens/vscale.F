*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VSCALE (B,ALPHA,A,N)
C
C CERN PROGLIB# F121    VSCALE          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),B(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= B(I) * ALPHA
C
      RETURN
      END
