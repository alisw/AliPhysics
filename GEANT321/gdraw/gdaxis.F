*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GDAXIS(X0,Y0,Z0,AXSIZ)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the axis system at a given                          *
C.    *       point (X0,Y0,Z0) with axis of length AXSIZ;              *
C.    *       by this way the user can easily realize                  *
C.    *       wich view is in effect in the picture                    *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDSPEC                      *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3,4),U(4),V(4),XL(3)
C
C.
C.    ------------------------------------------------------------------
C.
      DO 10 I=1,4
         X(1,I)=X0
         X(2,I)=Y0
         X(3,I)=Z0
   10 CONTINUE
C
      X(1,2)=X0+AXSIZ
      X(2,3)=Y0+AXSIZ
      X(3,4)=Z0+AXSIZ
      CALL GDFR3D(X,4,U,V)
C
      DO 20 I=1,3
         XL(I)=SQRT((U(I+1)-U(1))**2+(V(I+1)-V(1))**2)
   20 CONTINUE
C
      XLMAX=MAX(XL(1),XL(2))
      XLMAX=MAX(XLMAX,XL(3))
      ARSIZ=XLMAX/7.0
C
      CALL GDARR(U(1),V(1),U(2),V(2),ARSIZ,'<X')
      CALL GDARR(U(1),V(1),U(3),V(3),ARSIZ,'<Y')
      CALL GDARR(U(1),V(1),U(4),V(4),ARSIZ,'<Z')
C
      END
