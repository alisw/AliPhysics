*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFCOOR(X,IAXIS,C)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    ROUTINE TO COMPUTE THE COORDINATE C CORRESPONDING TO THE   *
C.    *    AXIS IAXIS FOR THE POINT X(1-3).                           *
C.    *        IAXIS = 1 IS THE USUAL X AXIS ------>  X(1)            *
C.    *        IAXIS = 2 IS THE USUAL Y AXIS ------>  X(2)            *
C.    *        IAXIS = 3 IS THE USUAL Z AXIS ------>  X(3)            *
C.    *        IAXIS = 4 IS THE X,Y PLANE R VALUE.                    *
C.    *        IAXIS = 5 IS THE FULL X,Y,Z R VALUE.                   *
C.    *        IAXIS = 6 IS AZIMUTH PHI (ANGLE IN X,Y PLANE).         *
C.    *        IAXIS = 7 IS THE POLAR ANGLE THETA (R,Z PLANE).        *
C.    *                                                               *
C.    *    ==>Called by : GFTRAC, GINVOL, GMEDIA, GNEXT, GTMEDI,      *
C.    *                  GTNEXT                                       *
C.    *         Author  A.McPherson  *********                        *
C.    *                                                               *
C.    *****************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
      DIMENSION X(*)
C.
C.         ----------------------------------------------------
C.
      IF(IAXIS.LE.0.OR.IAXIS.GE.8) THEN
         WRITE(CHMAIL,10000) IAXIS
         CALL GMAIL(0,0)
*
      ELSEIF(IAXIS.LE.3) THEN
         C=X(IAXIS)
*
      ELSEIF(IAXIS.EQ.4) THEN
         C=X(1)*X(1)+X(2)*X(2)
         IF(C.GT.0.0) C=SQRT(C)
*
      ELSEIF(IAXIS.EQ.5) THEN
         C=X(1)*X(1)+X(2)*X(2)+X(3)*X(3)
         IF(C.GT.0.0) C=SQRT(C)
*
      ELSEIF(IAXIS.EQ.6) THEN
         C=0.0
         IF(ABS(X(1)).GE.1.0E-10.AND.ABS(X(2)).GE.1.0E-10) THEN
            C=ATAN2(X(2),X(1))*RADDEG
         ELSEIF(ABS(X(4)).GE.1.0E-10.AND.ABS(X(5)).GE.1.0E-10) THEN
            C=ATAN2(X(5),X(4))*RADDEG
         ENDIF
         IF(C.LT.0.0) C=C+360.0
*
      ELSEIF(IAXIS.EQ.7) THEN
         C=90.0
         IF(ABS(X(3)).LT.1.0E-10) GO TO 999
         R=X(1)*X(1)+X(2)*X(2)
         IF(R.GT.0.0) R=SQRT(R)
         C=ATAN2(R,X(3))*RADDEG
         IF(C.LT.0.0) C=C+180.0
         IF(C.EQ.0.0.AND.X(3).LE.0.0) C=180.0
*
      ENDIF
10000 FORMAT(' ***** GFCOOR: IAXIS =',I5,' NOT SUPPORTED *****')
  999 END
