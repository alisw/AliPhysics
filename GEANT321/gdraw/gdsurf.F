*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDSURF(N)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the N-th surface described by the 4 space points    *
C.    *       stored in vector SURF starting from address 12*(N-1)+1   *
C.    *                                                                *
C.    *    ==>Called by : (not called yet)                             *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION SURF(2)
      EQUIVALENCE (WS(1),SURF(1))
      DIMENSION P(3,4),U(5),V(5)
C.
C.    ------------------------------------------------------------------
C.
      IADR=12*(N-1)+1
C
      DO 10  I=1,4
         P(1,I)=SURF(IADR)
         P(2,I)=SURF(IADR+1)
         P(3,I)=SURF(IADR+2)
         IADR=IADR+3
   10 CONTINUE
C
      CALL GDFR3D(P,4,U,V)
      U(5)=U(1)
      V(5)=V(1)
      CALL GDRAWV(U,V,5)
      END
