*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNOHYP (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH HYPErboloid      *
C.    *      VOLUME, FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6)  *
C.    *                                                                *
C.    *       PAR   (input)  : volume parameters                       *
C.    *                        inner radius                            *
C.    *                        outer radius                            *
C.    *                        half length in z                        *
C.    *                        straw stereo angle in degrees           *
C.    *                        r**2 = (z*tan(theta))**2 + a**2         *
C.    *                                                                *
C.    *       IACT  (input)  : action flag                             *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, and SNXT only if SNEXT .GT.new SAFE *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *       SNEXT (input)  : see IACT = 1                            *
C.    *       SNXT  (output) : distance to volume boundary             *
C.    *       SAFE  (output) : shortest distance to any boundary       *
C.    *                                                                *
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Authors  M.J.Corden, A.Palounek                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6),PAR(4),ENDR(2)
*
*     ------------------------------------------------------------------
*
      SNXT = BIG
      TANTHS= (TAN(PAR(4)*DEGRAD))**2
      R2 = X(1)*X(1)+X(2)*X(2)
      ENDR(1)=SQRT(PAR(1)**2+PAR(3)**2*TANTHS)
      ENDR(2)=SQRT(PAR(2)**2+PAR(3)**2*TANTHS)
      SAF3=ABS(X(3))-PAR(3)
*
      IF(IACT.EQ.3) GO TO 10
*
*     compute SAFE from escribed cylinders
*
      R  = SQRT(R2)
      SAF1=PAR(1)-R
      SAF2=R-ENDR(2)
      SAFE = MAX(0.,SAF1,SAF2,SAF3)
*
      IF(IACT.EQ.0) GO TO 999
      IF(IACT.EQ.1.AND.SAFE.GT.SNEXT) GO TO 999
*
   10 CONTINUE
*
* *** Compute SNXT
*
      IF((SAF3.GT.0.).AND.(X(3)*X(6).GE.0.)) GO TO 999
C
C       Compute intercept with inner & outer surfaces.
C
      A = X(4)**2 + X(5)**2 - X(6)**2*TANTHS
      B = X(1)*X(4) + X(2)*X(5) - X(3)*X(6)*TANTHS
C
      DO 30 NCYL =1,2
         C = R2 - X(3)**2*TANTHS - PAR(NCYL)**2
         SURD = B**2 - A*C
         IF(SURD.LE.0.0) GO TO 30
         SURD=SQRT(SURD)
*
         DO 20 IR=-1,+1,2
            S=(-B+IR*SURD)/A
            IF((S.LT.0.).OR.(S.GT.SNXT)) GO TO 20
*       Intersection point
            ZI = X(3) + S*X(6)
*       Check Z limits
            IF(ABS(ZI).GT.PAR(3)) GO TO 20
*       Intersection with hyperboloid ok
            SNXT = S
   20    CONTINUE
   30 CONTINUE
*
* *** Intersection with Z end planes
*
      IF(SAF3.LE.0.) GO TO 999
      S    = SAF3/ABS(X(6))
*     Intersection point
      XI  = X(1) + S*X(4)
      YI  = X(2) + S*X(5)
*     Check R limits
      RI  = SQRT(XI*XI + YI*YI)
      IF(RI.LT.ENDR(1).OR.RI.GT.ENDR(2)) GO TO 999
*     Intersection with Z end plane ok
      SNXT  = S
*
  999 END
