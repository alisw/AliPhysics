*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DFRDH2(F,G,X,WS,IDIM,N)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='DFRDH2')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RFRDH2(F,G,X,WS,IDIM,N)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RFRDH2')
#endif
C
      DIMENSION WS(IDIM,*)

      SUM=F(X)
      DO 2 J = 1,N
    2 SUM=SUM+WS(J,N+4)*G(X,WS(J,N+3))*WS(J,N+1)
      GO TO 9

#if defined(CERNLIB_DOUBLE)
      ENTRY DFRDH3(H,WS,IDIM,N)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY RFRDH3(H,WS,IDIM,N)
#endif

      SUM=0
      DO 3 J = 1,N
    3 SUM=SUM+WS(J,N+4)*H(WS(J,N+3))*WS(J,N+1)
#if defined(CERNLIB_DOUBLE)
    9 DFRDH2=SUM
#endif
#if !defined(CERNLIB_DOUBLE)
    9 RFRDH2=SUM
#endif
      RETURN
      END
