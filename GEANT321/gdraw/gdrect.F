*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRECT(P1,P2,P3,P4)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a rectangle of points P1,P2,P3,P4                   *
C.    *                                                                *
C.    *    ==>Called by : GDRAWS                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
CM
#include "geant321/gcdraw.inc"
CM
      DIMENSION P1(3),P2(3),P3(3),P4(3),U(5),V(5),X(3,4)
      DIMENSION Q1(3),Q2(3),Q3(3),Q4(3)
      EQUIVALENCE (Q1(1),X(1,1)),(Q2(1),X(1,2))
      EQUIVALENCE (Q3(1),X(1,3)),(Q4(1),X(1,4))
C.
C.    ------------------------------------------------------------------
C.
      DO 10 I=1,3
         Q1(I)=P1(I)
         Q2(I)=P2(I)
         Q3(I)=P3(I)
         Q4(I)=P4(I)
   10 CONTINUE
C
      CALL GDFR3D(X,4,U,V)
C
CM             Short circuit in case of GMR
C
      IF ( IGMR .NE. 0 ) GO TO 999
CM
      U(5)=U(1)
      V(5)=V(1)
      CALL GDRAWV(U,V,5)
C
  999 RETURN
      END
