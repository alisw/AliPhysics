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
      SUBROUTINE GNCTUB (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH 'CTUB'           *
C.    *      VOLUME FROM INSIDE POINT X(1-3) ALONG DIRECTION X(4-6)    *
C.    *                                                                *
C.    *       PAR   (input)  : volume parameters                       *
C.    *       IACT  (input)  : action flag                             *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, and SNXT only if SNEXT .GT.new SAFE *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *       SNEXT (input)  : see IACT = 1                            *
C.    *       SNXT  (output) : distance to volume boundary             *
C.    *       SAFE  (output) : shortest distance to any boundary       *
C.    *                                                                *
C.    *    ==>Called by : GNEXT, GNPCON, GTNEXT                        *
C.    *         Authors A.McPherson    ********                        *
C.    *       MODIFICATION LOG.                                        *
C.    *       18-July-89 M.Guckes modifications due to GEANG 3.13      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      REAL X(6),PAR(11)
C
C-------------------------------------------------------------
C
      SNXT = BIG
      R2   = X(1)*X(1) + X(2)*X(2)
      R    = SQRT (R2)
      IF(PAR(5).GE.360.) THEN
         IFLAG = 1
      ELSE
         IFLAG = 2
      ENDIF
*
      SAFZ1  = (-PAR(3)-X(3) )*PAR(8)-X(1)*PAR(6)-X(2)*PAR(7)
      SAFZ2  = (PAR(3)-X(3) )*PAR(11)-X(1)*PAR(9)-X(2)*PAR(10)
      IF (PAR(1).NE.0.) THEN
         SAFR1  = R - PAR(1)
      ELSE
         SAFR1  = BIG
      ENDIF
      SAFR2  = PAR(2) - R
*
      IF (IACT .LT. 3) THEN
*
* *** Compute safety-distance 'SAFE' (P.Weidhaas)
*
         SAFSEG = BIG
         IF (IFLAG .EQ. 2) THEN
*
*     In addition to the radial distances (SAFR1 and SAFR2) and the
*     axial distances (SAFZ1 and SAFZ2) we compute here the distance
*     to the PHI-segment boundary that is closest to the point:
*
*     For each PHI-boundary we find the distance from the given
*     point to the outer (at RMAX) point of the segment boundary
*     (DISTS1 and DISTS2, resp.). If DISTS1 < DISTS2, we define
*     SAFSEG to be the distance to segment PHI1. Else we set
*     SAFSEG to be the distance to segment PHI2.
*
            PHI1 = PAR(4) * DEGRAD
            PHI2 = PAR(5) * DEGRAD
*
            COSPH1 = COS (PHI1)
            COSPH2 = COS (PHI2)
            SINPH1 = SIN (PHI1)
            SINPH2 = SIN (PHI2)
*
* ***   Get coordinates of outer endpoints (at ROUT) of both segments.
*
            XS1 = PAR(2) * COSPH1
            YS1 = PAR(2) * SINPH1
            XS2 = PAR(2) * COSPH2
            YS2 = PAR(2) * SINPH2
*
* ***   Get distances (squared) from the given point to each endpoint.
*
            DISTS1 = (X(1) - XS1)**2 + (X(2) - YS1)**2
            DISTS2 = (X(1) - XS2)**2 + (X(2) - YS2)**2
*
* ***   Get distance to that PHI-segment whose endpoint
* ***   is closest to the given point.
*
            IF (DISTS1 .LE. DISTS2) THEN
               SAFSEG = ABS(X(1) * SINPH1 - X(2) * COSPH1)
            ELSE
               SAFSEG = ABS(X(1) * SINPH2 - X(2) * COSPH2)
            ENDIF
         ENDIF
*
         SAFE = MIN (SAFZ1, SAFZ2, SAFR1, SAFR2, SAFSEG)
         IF (IACT .EQ. 0) GO TO 999
         IF (IACT .EQ. 1 .AND. SNEXT .LT. SAFE) GO TO 999
      ENDIF
*
* ***   Compute intersection with z-boundaries
*
      V1 = X(4)*PAR(6)+X(5)*PAR(7)+X(6)*PAR(8)
      SZ1 = SAFZ1/V1
      V2 = X(4)*PAR(9)+X(5)*PAR(10)+X(6)*PAR(11)
      SZ2 = SAFZ2/V2
      IF( SZ1 .GT. 0. ) THEN
         SNXT = SZ1
      ELSE
         SNXT = BIG
      ENDIF
      IF( SZ2 .GT. 0.0 .AND. SZ2 .LT. SNXT) SNXT = SZ2
      SZ=SNXT
*
      IF (ABS(X(6)).LT.1.)THEN
*
* ***      Compute z-intercept with inner cylinder.
*
         BA2=-1.
         IF(PAR(1).GT.0.)THEN
            RMIN2  = PAR(1)*PAR(1)
            ZP2    = 1./(1.-X(6)*X(6))
            BA     = (X(4)*X(1)+X(5)*X(2))*ZP2
            BA2    = BA*BA
            CA     = (R2-RMIN2)*ZP2
            DIS2   = BA2-CA
            IF (DIS2.GE.0.)THEN
               XSIN    = -BA-SQRT(DIS2)
               IF (XSIN.GE.0.)THEN
                  IF(XSIN.LT.SNXT)SNXT = XSIN
                  GO TO 10
               ENDIF
            ENDIF
         ENDIF
*
*  ***    Compute z-intercept with outer cylinder.
*
         RMAX2  = PAR(2)*PAR(2)
         XZ     = X(1) + X(4)*SZ
         YZ     = X(2) + X(5)*SZ
         IF (XZ*XZ+YZ*YZ.GT.RMAX2)THEN
            IF(BA2.LT.0.)THEN
               ZP2    = 1./(1.-X(6)*X(6))
               BA     = (X(4)*X(1)+X(5)*X(2))*ZP2
               BA2    = BA*BA
            ENDIF
            CA     = (R2-RMAX2)*ZP2
            DIS2   = BA2-CA
            IF (DIS2.GE.0.)THEN
               SRMAX = -BA+SQRT(DIS2)
               IF(SRMAX.LT.SNXT)SNXT=SRMAX
            ENDIF
         ENDIF
      ENDIF
*
   10 IF(IFLAG.EQ.2) THEN
*
*     =======>PHI segmented tube
*             We have checked the radius and Z.
*             now check PHI.
*
         DPSGN=X(1)*X(5)-X(2)*X(4)
*
*             Tells us which way its going.
*
         PHI2=PAR(5)
         IF(DPSGN.LT.0.0) PHI2=PAR(4)
         PH2=PHI2*DEGRAD
*
*             Have set up the limit.
*
         TSGN=SIN(PH2)
         TCSG=COS(PH2)
         DX45=X(4)*TSGN-X(5)*TCSG
         IF(DX45.EQ.0.)GO TO 999
         SN1=(X(2)*TCSG-X(1)*TSGN)/DX45
*
*             Distance until tangents are right.
*
         IF(SN1.LT.0.0) GO TO 999
         IF(ABS(TSGN).GT.1.E-6.AND.(X(2)+SN1*X(5))*TSGN .LT.0.)GO TO
     +   999
         IF(ABS(TCSG).GT.1.E-6.AND.(X(1)+SN1*X(4))*TCSG .LT.0.)GO TO
     +   999
*
*             Have checked that the distance is +VE and that the
*             SINE is the right sign.
*
         IF(SN1.LT.SNXT) SNXT=SN1
      END IF
*
  999 END
