*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE RSRTNT(K,N,A,B,C,U1,V1,RES,LRL)
#include "gen/def64.inc"
     +  DRS,D
      LOGICAL LRL

      SROUND(D)=D+(D-SNGL(D))
      CALL DSRTNT(K,N,DBLE(A),DBLE(B),DBLE(C),DBLE(U1),DBLE(V1),DRS,LRL)
      RES=SROUND(DRS)
      RETURN
      END
#endif
