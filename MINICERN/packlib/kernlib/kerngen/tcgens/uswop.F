*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE USWOP (A,B,N)
C
C CERN PROGLIB# V303    USWOP           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71 RKB
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*),B(*), X
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*),B(*)
#endif
C
      IF (N.EQ.0) RETURN
         DO 10 I=1,N
      X=A(I)
      A(I)=B(I)
   10 B(I)=X
      RETURN
      END
