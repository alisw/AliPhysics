*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDLINE(P1,P2)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a line between points P1 and P2                     *
C.    *                                                                *
C.    *    ==>Called by : GDLCYL, GDRAWS                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
CM
#include "geant321/gcdraw.inc"
CM
      DIMENSION P1(3),P2(3),U(2),V(2),X(3,2)
      DIMENSION Q1(3),Q2(3)
      EQUIVALENCE (Q1(1),X(1,1)),(Q2(1),X(1,2))
C.
C.    ------------------------------------------------------------------
C.
      DO 10 I=1,3
         Q1(I)=P1(I)
         Q2(I)=P2(I)
   10 CONTINUE
C
      CALL GDFR3D(X,2,U,V)
C
CM             Short circuit in case of GMR
C
      IF ( IGMR .EQ. 0 ) CALL GDRAWV(U,V,2)
CM
      END
