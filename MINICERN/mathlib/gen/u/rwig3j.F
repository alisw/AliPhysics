*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RWIG3J(A,B,C,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (D)
 
      SROUND(D)=D+(D-SNGL(D))
      RWIG3J=SROUND(
     1 DWIG3J(DBLE(A),DBLE(B),DBLE(C),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
 
      ENTRY RCLEBG(A,B,C,X,Y,Z)
      RCLEBG=SROUND(
     1 DCLEBG(DBLE(A),DBLE(B),DBLE(C),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
 
      ENTRY RWIG6J(A,B,C,X,Y,Z)
      RWIG6J=SROUND(
     1 DWIG6J(DBLE(A),DBLE(B),DBLE(C),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
 
      ENTRY RRACAW(A,B,C,X,Y,Z)
      RRACAW=SROUND(
     1 DRACAW(DBLE(A),DBLE(B),DBLE(C),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
 
      ENTRY RJAHNU(A,B,C,X,Y,Z)
      RJAHNU=SROUND(
     1 DJAHNU(DBLE(A),DBLE(B),DBLE(C),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
 
      ENTRY RWIG9J(A,B,C,P,Q,R,X,Y,Z)
      RWIG9J=SROUND(DWIG9J(DBLE(A),DBLE(B),DBLE(C),DBLE(P),DBLE(Q),
     1                     DBLE(R),DBLE(X),DBLE(Y),DBLE(Z)))
      RETURN
      END

#endif
