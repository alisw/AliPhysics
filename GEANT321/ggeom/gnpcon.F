*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNPCON (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'PCON' VOLUME,  *
C.    *       FROM INSIDE POINT X(1-3) ALONG DIRECTION X(4-6)          *
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
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Author  A.McPherson,  P.Weidhaas  *********            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6), PAR(9),  XT(6), PT(7)
 
      EQUIVALENCE (PT(1), DZ),  (PT(2), PT2), (PT(3), PT3)
      EQUIVALENCE (PT(4), PT4), (PT(5), PT5), (PT(6), PT6)
      EQUIVALENCE (PT(7), PT7)
 
      EQUIVALENCE (XT(3), XT3)
C.
C.   ---------------------------------------------------
C.
 
      SNXT = BIG
      R2  = X(1)*X(1) + X(2)*X(2)
      R   = SQRT (R2)
 
      NZ  = PAR(3)
      ZMIN  = PAR(4)
      ZMAX  = PAR(3*NZ+1)
      SAFZ1 = X(3) - ZMIN
      SAFZ2 = ZMAX - X(3)
      SAFEZ = MIN (SAFZ1, SAFZ2)
 
C
C......  First determine in which z-segment the particle is located.
C
      DO 10 JPH=7, 3*(NZ-1)+1,  3
        IF (X(3) .LE. PAR(JPH)) THEN
           IPH=JPH
           GO TO 20
        ENDIF
   10 CONTINUE
      IPH=3*NZ+1
 
   20 CONTINUE
C
C......  The particle is in the segment bounded by z-planes at
C......  Z1=PAR(IPL) and Z2=PAR(IPH), i.e.,  Z1 < X(3) < Z2.
C
C......  Set parameters for this segment and translate z-coordinate
C......  of point relative to center of this segment.This is done in
C......  preparation of invoking the algorithms used in "GNTUBE" and
C......  "GNCONE" (which for reasons of efficiency and clarity are
C......  implemented inline).
C
      IPL  = IPH - 3
      DZ    = 0.5 * (PAR(IPH) - PAR(IPL))
      PT2 = PAR(IPL+1)
      PT3 = PAR(IPL+2)
      PT4 = PAR(IPH+1)
      PT5 = PAR(IPH+2)
      PT6 = PAR(1)
      PT7 = PAR(1) + PAR(2)
      IF (PT7 .GT. 360.0)  PT7  = PT7 - 360.0
 
      XT3 = X(3) - 0.5 * (PAR(IPL) + PAR(IPH))
 
      XT(1)  = X(1)
      XT(2)  = X(2)
      XT(4)  = X(4)
      XT(5)  = X(5)
      XT(6)  = X(6)
 
      IND = 2
      IF (PAR(2) .EQ. 360.0)  IND = 1
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
        SAFZ1   = DZ - ABS(XT3)
        SAFEZ   = MIN (SAFEZ,SAFZ1)
        SAFSEG  = BIG
C
C......  Determine whether the segment is a tube or a cone.
C
 
        IF (PT2 .NE. PT4) GO TO 50
        IF (PT3 .NE. PT5) GO TO 50
 
C*********************************************************
C
C......  The segment is a tube: invoke the algorithm
C......  from routine "GNTUBE" inline to get "SAFER".
C
C*********************************************************
 
        SAFR1  = R - PT2
        SAFR2  = PT3 - R
        SAFER  = MIN (SAFR1, SAFR2)
 
        IF (IND .EQ. 2) GO TO 70
 
        GO TO 100
 
 
   50   CONTINUE
 
C*********************************************************
C
C......  The segment is a cone: invoke the algorithm
C......  from routine "GNCONE" inline to get "SAFER".
C
C*********************************************************
 
        SAFZ2  = DZ + XT3
        ZLENI  = 0.5 / DZ
 
C......  Compute radial distance to inner wall.
 
        FACT  = (PT4 - PT2) * ZLENI
        R1  = PT2 + FACT * SAFZ2
        SAFR1 = (R - R1) / SQRT(1.0 + FACT*FACT)
 
C......  Compute radial distance to outer wall.
 
        FACT  = (PT5 - PT3) * ZLENI
        R2  = PT3 + FACT * SAFZ2
        SAFR2 = (R2 - R) / SQRT(1.0 + FACT*FACT)
        SAFER = MIN (SAFR1, SAFR2)
 
        IF (IND .EQ. 1) GO TO 100
 
   70   CONTINUE
 
C********************************************************************
C......  Here we handle the case of a PHI-segment of a tube or cone.
C......  In addition to the radial distances (SAFR1, SAFR2) and the
C......  axial distances (SAFZ1, SAFZ2) we compute here the distance
C......  to the PHI-segment boundary that is closest to the point.
C
C......  For each PHI-boundary we find the distance from the given
C......  point to the outer (at R2) point of the segment boundary
C......  (DISTS1 and DISTS2, resp.). If DISTS1 < DISTS2, we define
C......  "SAFSEG" to be the distance to segment PHI1, else we set
C......  "SAFSEG" to be the distance to segment PHI2.
C*********************************************************************
 
          PHI1  = PT6 * DEGRAD
          PHI2  = PT7 * DEGRAD
          IF (PHI2 .LT. PHI1)  PHI2 = PHI2 + TWOPI
 
          COSPH1  = COS (PHI1)
          COSPH2  = COS (PHI2)
          SINPH1  = SIN (PHI1)
          SINPH2  = SIN (PHI2)
 
C......  Get coordinates of outer endpoints (at R2) of both PHI-segments.
 
          XS1  = R2 * COSPH1
          YS1  = R2 * SINPH1
          XS2  = R2 * COSPH2
          YS2  = R2 * SINPH2
 
C......  Get distances (squared) from given point to each endpoint.
 
          DISTS1 = (X(1) - XS1)**2  +  (X(2) - YS1)**2
          DISTS2 = (X(1) - XS2)**2  +  (X(2) - YS2)**2
 
C......  Get distance to that PHI-segment whose endpoint
C......  is closest to the given point.
 
          IF (DISTS1 .LE. DISTS2) THEN
            SAFSEG = ABS(X(1) * SINPH1 - X(2) * COSPH1)
          ELSE
            SAFSEG = ABS(X(1) * SINPH2 - X(2) * COSPH2)
          ENDIF
 
 
  100   CONTINUE
 
        SAFE  = MIN (SAFEZ, SAFER, SAFSEG)
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
*
*     Avoid rounding effects induced by translation ********************
*
      IF (ABS(XT(3)).GE.DZ) XT(3) = (1.-0.000001)*XT(3)
*
      IF (PT2 .NE. PT4) GO TO 200
      IF (PT3 .NE. PT5) GO TO 200
 
      DELZ  = DZ
      PT(1) = PT2
      PT(2) = PT3
      PT(3) = DELZ
      PT(4) = PT6
      PT(5) = PT7
 
      CALL GNTUBE (XT, PT, 3, IND, SNEXT, SNXT, TSAFE)
      GO TO 999
 
  200 CONTINUE
      CALL GNCONE (XT, PT, 3, IND, SNEXT, SNXT, TSAFE)
 
  999 CONTINUE
      END
 
