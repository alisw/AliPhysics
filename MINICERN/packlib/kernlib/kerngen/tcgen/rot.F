*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE ROT (A,TH,B)
C
C CERN PROGLIB# F118    ROT             .VERSION KERNFOR  1.0   661115
C ORIG. 15/11/66
C
      COMMON /SLATE/B1,B2,COSTH,SINTH,XQX(36)
      DIMENSION A(3),B(3)
C
      COSTH= COS(TH)
      SINTH= SIN(TH)
      B1  = A(1)*COSTH - A(2)*SINTH
      B2  = A(1)*SINTH + A(2)*COSTH
      B(3)= A(3)
      B(2)= B2
      B(1)= B1
      RETURN
      END
