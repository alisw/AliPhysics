*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:27  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE PROC1(A,QR,M,N,IM)
      DIMENSION A(IM,N),QR(M,N)
C
      DO 5 I=1,M
      DO 5 J=1,N
 5    QR(I,J)=A(I,J)
      RETURN
      END
