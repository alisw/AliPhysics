*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UCOPIV (A,B,N)
C
C CERN PROGLIB# V301    UCOPIV          .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*), B(*)
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*), B(*)
#endif
C
      IF (N.EQ.0) RETURN
      I = N
      DO 12 J=1,N
      B(I)= A(J)
   12 I = I - 1
      RETURN
      END
