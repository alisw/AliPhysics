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
      SUBROUTINE GNHYPE (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH  HYPErboloid     *
C.    *      VOLUME FROM INSIDE POINT X(1-3) ALONG DIRECTION X(4-6)    *
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
C.    *         Authors M.J. Corden,  A.Palounek                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6),PAR(4)
 
C-------------------------------------------------------------
 
      SNXT = BIG
      R2   = X(1)*X(1) + X(2)*X(2)
      TANTHS=(TAN(PAR(4)*DEGRAD))**2
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
      SAFZ1  = PAR(3) - X(3)
      SAFZ2  = PAR(3) + X(3)
 
      IF (IACT .LT. 3) THEN
 
C    tube version:
C      R    = SQRT (R2)
C      IF (PAR(1).NE.0.) THEN
C         SAFR1  = R - PAR(1)
C      ELSE
C         SAFR1  = BIG
C      ENDIF
C      SAFR2  = PAR(2) - R
C
C        SAFE  = MIN (SAFZ1, SAFZ2, SAFR1, SAFR2)
C
C    simple, safe choice for hyperboloid for now:
      SAFE=0.
C
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
C
C             Compute intersection with z-boundaries.
C
      IF(X(6).GT.0.)THEN
         SNXT=SAFZ1/X(6)
      ELSEIF(X(6).LT.0.)THEN
         SNXT=-SAFZ2/X(6)
      ENDIF
C
C     Compute intercepts with inner & outer curved surfaces:
C
      A = X(4)**2 + X(5)**2 - X(6)**2*TANTHS
      B = X(1)*X(4) + X(2)*X(5) - X(3)*X(6)*TANTHS
 
      DO 38 NCYL =1,2
        C = R2 - X(3)**2*TANTHS - PAR(NCYL)**2
        SURD = B**2 - A*C
        IF(SURD.GT.0.) THEN
          SURD=SQRT(SURD)
          DIST1=(-B+SURD)/A
          IF(DIST1.GT.0.) SNXT=MIN(SNXT,DIST1)
          DIST2=(-B-SURD)/A
          IF(DIST2.GT.0.) SNXT=MIN(SNXT,DIST2)
C       not necessary to compute intercept with outer surface
C       if there is a positive intercept with inner surface:
          IF(DIST1.GT.0. .OR. DIST2.GT.0.) GO TO 999
        ENDIF
  38  CONTINUE
C
 999  CONTINUE
      END
