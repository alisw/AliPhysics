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
      SUBROUTINE GITRAN(X,DX,IROT,XNEW)
C.
C.    ******************************************************************
C     *                                                                *
C     *        ROUTINE TO APPLY THE TRANSFORMATION GIVEN BY            *
C     *        THE POINT DX AND THE ROTATION MATRIX IN THE GEANT       *
C     *        ROTATION MATRIX BANKS WITH INDEX NUMBER IROT TO THE     *
C     *        POINT X, STORING THE RESULT AS THE POINT XNEW.          *
C     *                                                                *
C.    *    ==>Called by : GFTRAC,GINVOL,GMEPOS,GNEXT,GTMEDI,GTNEXT     *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION X(3),DX(3),XNEW(3)
C.
C.    ------------------------------------------------------------------
C.
      IF(IROT.EQ.0) THEN
         XNEW(1) = X(1) - DX(1)
         XNEW(2) = X(2) - DX(2)
         XNEW(3) = X(3) - DX(3)
*
      ELSE
         XL1 = X(1) - DX(1)
         XL2 = X(2) - DX(2)
         XL3 = X(3) - DX(3)
         JR=LQ(JROTM-IROT)
         XNEW(1) = XL1*Q(JR+1) + XL2*Q(JR+2) + XL3*Q(JR+3)
         XNEW(2) = XL1*Q(JR+4) + XL2*Q(JR+5) + XL3*Q(JR+6)
         XNEW(3) = XL1*Q(JR+7) + XL2*Q(JR+8) + XL3*Q(JR+9)
*
      ENDIF
      END
