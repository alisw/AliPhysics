*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCIRC(R,Z0)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a circle of radius R centered on Z-axis,            *
C.    *       at a distance of Z0 from the origin,                     *
C.    *       on a plane parallel to X-Y plane                         *
C.    *                                                                *
C.    *    ==>Called by : GDRAWS                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
      DIMENSION X(3,41),U(41),V(41)
C.
C.    ------------------------------------------------------------------
C.
      IF (R.LT.0.0001) GO TO 999
C
      DO 10 I=1,40
         X(1,I)=R*GCOS(I)
         X(2,I)=R*GSIN(I)
         X(3,I)=Z0
   10 CONTINUE
C
      CALL GDFR3D(X,40,U,V)
C
CM             Short circuit in case of GMR
C
      IF ( IGMR .NE. 0 ) GO TO 999
CM
      U(41)=U(1)
      V(41)=V(1)
      CALL GDRAWV(U,V,41)
C
  999 RETURN
      END
