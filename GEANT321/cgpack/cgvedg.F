*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGVEDG(NT,P1,P2,IVIS)
************************************************************************
*                                                                      *
*     Name: CGVEDG                                                     *
*     Author: E. Chernyaev                       Date:    08.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Compute ALL about edge:                                *
*               parametric equations X=X1+T(X2-X1), Y=Y1+T(Y2-Y1)      *
*               edge line equation   Ax+By+C=0                         *
*               edge min-max                                           *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: NT - number of transformation to screen coordinates       *
*            P1(3) - 1-st edge point                                   *
*            P2(3) - 2-nd edge point                                   *
*                                                                      *
*     Output: IVIS = -1 for very small edges                           *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgdelt.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
      REAL    P1(3),P2(3)
*-
      DO 100 I=1,3
        AA(I) = P1(1)*TSCRN(1,I,NT) + P1(2)*TSCRN(2,I,NT) +
     +          P1(3)*TSCRN(3,I,NT) + TSCRN(4,I,NT)
        BB(I) = P2(1)*TSCRN(1,I,NT) + P2(2)*TSCRN(2,I,NT) +
     +          P2(3)*TSCRN(3,I,NT) + TSCRN(4,I,NT)
        DELTA(I) = BB(I) - AA(I)
  100   CONTINUE
*
      IVIS   =-1
      AXDELT=ABS(XDELT)
      AYDELT=ABS(YDELT)
      IF ((AXDELT.LE.EESCR).AND.(AYDELT.LE.EESCR))GOTO 999
      IVIS   = 0
*           C O M P U T E   E D G E   E Q U A T I O N
      Q        = XDELT*XDELT + YDELT*YDELT
      ABC(1,1) = XDELT/Q
      ABC(2,1) = YDELT/Q
      ABC(3,1) =-(XDELT*XA + YDELT*YA)/Q
      XXXQ        = SQRT(Q)
      ABC(1,2) =-YDELT/XXXQ
      ABC(2,2) = XDELT/XXXQ
      ABC(3,2) = (YDELT*XA - XDELT*YA)/XXXQ
*           C O M P U T E   E D G E   M I N - M A X
      DEDGE(1) =-MAX(XA,XB)
      DEDGE(2) =-MAX(YA,YB)
      DEDGE(3) = MIN(XA,XB)
      DEDGE(4) = MIN(YA,YB)
      DEDGE(5) = MIN(ZA,ZB)
*
      NPART  = 0
      TDEL   = EESCR / XXXQ
  999 RETURN
      END
