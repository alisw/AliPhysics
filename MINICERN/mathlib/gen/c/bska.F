*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE BSKA(X,IA,JA,NL,B)
      DOUBLE PRECISION D,DB(0:100)
      DIMENSION B(0:*)
 
      SROUND(D)=D+(D-SNGL(D))
 
      CALL DBSKA(DBLE(X),IA,JA,NL,DB)
      GO TO 9
 
      ENTRY EBSKA(X,IA,JA,NL,B)
 
      CALL DEBKA(DBLE(X),IA,JA,NL,DB)
    9 DO 1 N = 1,NL
    1 B(N)=SROUND(DB(N))
      RETURN
      END
#endif
