*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINROT(X,RMAT,XNEW)
C.
C.    ******************************************************************
C     *                                                                *
C     *        ROUTINE TO ROTATE THE VECTOR X BY THE INVERSE OF THE    *
C     *       ROTATION MATRIX RMAT AND STORE THE RESULT AS THE         *
C     *        VECTOR XNEW                                             *
C.    *                                                                *
C.    *    ==>Called by : GDFR3D,GDSSUR,GD3D3D,GFLCAR,GFLPHI,GFLRAD    *
C.    *                   GFLTHE                                       *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),XNEW(3),RMAT(10)
C.
C.    ------------------------------------------------------------------
C.
      IF(RMAT(10).NE.0.0) THEN
         XNEW(1) = X(1)*RMAT(1) + X(2)*RMAT(4) + X(3)*RMAT(7)
         XNEW(2) = X(1)*RMAT(2) + X(2)*RMAT(5) + X(3)*RMAT(8)
         XNEW(3) = X(1)*RMAT(3) + X(2)*RMAT(6) + X(3)*RMAT(9)
*
      ELSE
         XNEW(1) = X(1)
         XNEW(2) = X(2)
         XNEW(3) = X(3)
*
      ENDIF
      END
