*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GROT(X,RMAT,XNEW)
C.
C.    ******************************************************************
C.    *                                                                *
C     *            ROUTINE TO ROTATE POINT X BY MATRIX RMAT TO         *
C     *            FORM POINT XNEW.                                    *
C     *                                                                *
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),XNEW(3),RMAT(10)
C
      IF(RMAT(10).NE.0.0) THEN
         XNEW(1) = X(1)*RMAT(1) + X(2)*RMAT(2) + X(3)*RMAT(3)
         XNEW(2) = X(1)*RMAT(4) + X(2)*RMAT(5) + X(3)*RMAT(6)
         XNEW(3) = X(1)*RMAT(7) + X(2)*RMAT(8) + X(3)*RMAT(9)
*
      ELSE
         XNEW(1) = X(1)
         XNEW(2) = X(2)
         XNEW(3) = X(3)
*
      ENDIF
      END
