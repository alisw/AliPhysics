*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE VMATR (A,G,X,N,K)
C
C CERN PROGLIB# F121    VMATR           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),X(*),G(K,N)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,K
      X(I)= 0.
      DO 9 J= 1,N
    9 X(I)= X(I) + G(I,J)*A(J)
      RETURN
      END
