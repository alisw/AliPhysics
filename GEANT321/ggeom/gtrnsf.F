*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GTRNSF(X,DX,RMAT,XNEW)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       ROUTINE TO APPLY THE TRANSFORMATION GIVEN BY             *
C     *       THE POINT DX AND THE ROTATION MATRIX RMAT TO THE         *
C     *       POINT X, PLACING THE RESULT IN THE POINT XNEW.           *
C.    *                                                                *
C.    *    ==>Called by : GDFR3D, GFTRAC, GINVOL,GMEDIA, GNEXT, GSTRAC,*
C.    *                   GTMEDI, GTNEXT                               *
C.    *         Author  A.McPherson  *********                         *
C     *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),DX(3),RMAT(10),XNEW(3)
C
      IF(RMAT(10).EQ.0.)THEN
         XNEW(1) = X(1) - DX(1)
         XNEW(2) = X(2) - DX(2)
         XNEW(3) = X(3) - DX(3)
*
      ELSE
         XL1 = X(1) - DX(1)
         XL2 = X(2) - DX(2)
         XL3 = X(3) - DX(3)
         XNEW(1) = XL1*RMAT(1) + XL2*RMAT(2) + XL3*RMAT(3)
         XNEW(2) = XL1*RMAT(4) + XL2*RMAT(5) + XL3*RMAT(6)
         XNEW(3) = XL1*RMAT(7) + XL2*RMAT(8) + XL3*RMAT(9)
 
      ENDIF
      END
