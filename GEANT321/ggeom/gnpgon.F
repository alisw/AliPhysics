*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNPGON (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'PGON' VOLUME,  *
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "geant321/gconsp.inc"
 
      REAL X(6), PAR(50), SNEXT, SNXT, SAFE
      DIMENSION TPAR(50)
 
C-----------------------------------------------------------------
 
      SNXT = BIG
      SAFSEG = BIG
      R2 = X(1)*X(1) + X(2)*X(2)
      R  = SQRT (R2)
 
      PDIV  = PAR(2) / PAR(3)
      DELPHI  = PDIV * DEGRAD
      DPHI2 = 0.5 * DELPHI
      CSDPH2 = COS (DPHI2)
 
      TPAR(1)  = PAR(1)
      TPAR(2)  = PAR(2)
      TPAR(3)  = PAR(3)
      TPAR(4)  = PAR(4)
      NZ  = PAR(4)
 
      DO 5 I=1, NZ
        I3  = 3*(I-1)
        TPAR(5+I3) = PAR(5+I3)
        TPAR(6+I3) = PAR(6+I3) / CSDPH2
        TPAR(7+I3) = PAR(7+I3)
    5 CONTINUE
 
C**************************************************************************
C
C......  Here we start the logic from "GNPCON" (which for reasons of
C......  efficiency and clarity has been implemented inline).
C
C**************************************************************************
 
      ZMIN  = TPAR(5)
      ZMAX  = TPAR(3*NZ+2)
      SAFZ1 = X(3) - ZMIN
      SAFZ2 = ZMAX - X(3)
      SAFEZ = MIN (SAFZ1, SAFZ2)
 
C
C......  First determine in which segment the particle is located.
C
      DO 10 JPH=8, 3*NZ-1,  3
        IF (X(3) .LT. TPAR(JPH)) THEN
           IPH=JPH
           GO TO 20
        ENDIF
   10 CONTINUE
      IPH = 3*NZ+2
 
   20 CONTINUE
C
C......  The particle is in the segment bounded by z-planes at
C......  Z1=PAR(IPL) and Z2=PAR(IPH), i.e.,  Z1 < X(3) < Z2.
C
C......  Set parameters for this segment and translate z-coordinate
C......  of point relative to center of this segment. this is done in
C......  preparation of invoking the algorithms used in "GNTUBE" and
C......  "GNCONE" (which for reasons of efficiency and clarity are
C......  implemented inline).
C
      IPL  = IPH - 3
      DZ    = 0.5 * (TPAR(IPH) - TPAR(IPL))
      PT2 = TPAR(IPL+1)
      PT3 = TPAR(IPL+2)
      PT4 = TPAR(IPH+1)
      PT5 = TPAR(IPH+2)
      PT6 = TPAR(1)
      PT7 = TPAR(1) + TPAR(2)
      IF (PT7 .GT. 360.0)  PT7  = PT7 - 360.0
 
      XT3 = X(3) - 0.5 * (TPAR(IPL) + TPAR(IPH))
 
      SAFZ2  = DZ + XT3
      ZLENI  = 0.5 / DZ
 
      PHI1  = PT6 * DEGRAD
      PHI2  = PT7 * DEGRAD
      IF (PHI2.LE.PHI1) PHI2=PHI2+TWOPI
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
        IFLAG = 2
        IF (TPAR(2) .EQ. 360.0)  IFLAG = 1
        SAFZ1 = DZ - ABS(XT3)
        SAFEZ = MIN (SAFEZ,SAFZ1)
C
C......  Next determine whether the segment is a tube or a cone.
C
        IF (PT2 .NE. PT4) GO TO 50
        IF (PT3 .NE. PT5) GO TO 50
 
C*********************************************************
C
C......  The segment is a tube; invoke the algorithm
C......  from routine "GNTUBE" inline to get "SAFER".
C
C*********************************************************
 
        SAFR1  = BIG
        IF(PT2.GT.0.) SAFR1  = R - PT2
        SAFR2  = PT3 - R
        SAFER  = MIN (SAFR1, SAFR2)
 
        IF (IFLAG .EQ. 2) GO TO 70
 
        GO TO 100
 
 
   50   CONTINUE
 
C*********************************************************
C
C......  The segment is a cone; invoke the algorithm
C......  from routine "GNCONE" inline to get "SAFER".
C
C*********************************************************
 
 
C......  Compute radial distance to inner wall.
 
        IF(PT2+PT4.GT.0.) THEN
        FACT  = (PT4 - PT2) * ZLENI
        RAD1  = PT2 + FACT * SAFZ2
        SAFR1 = (R - RAD1) / SQRT(1.0 + FACT*FACT)
        ELSE
        SAFR1 = BIG
        ENDIF
 
C......  Compute radial distance to outer wall.
 
        FACT  = (PT5 - PT3) * ZLENI
        RAD2  = PT3 + FACT * SAFZ2
        SAFR2 = (RAD2 - R) / SQRT(1.0 + FACT*FACT)
 
        SAFER  = MIN (SAFR1, SAFR2)
 
        IF (IFLAG .EQ. 1) GO TO 100
 
   70   CONTINUE
 
C********************************************************************
C......  Here we handle the case of a phi-segment of a tube or cone.
C......  in addition to the radial distances (SAFR1, SAFR2) and the
C......  axial distances (SAFZ1, SAFZ2) we compute here the distance
C......  to the phi-segment boundary that is closest to the point:
C
C......  For each phi-boundary we find the distance from the given
C......  point to the outer (at R2) point of the segment boundary
C......  (DISTS1 and DISTS2, resp.). If DISTS1 < DISTS2, we define
C......  "SAFSEG" to be the distance to segment PHI1, else we set
C......  "SAFSEG" to be the distance to segment PHI2.
C*********************************************************************
 
 
          COSPH1  = COS (PHI1)
          SINPH1  = SIN (PHI1)
          COSPH2  = COS (PHI2)
          SINPH2  = SIN (PHI2)
 
C......  Get coordinates of outer endpoints (at R2) of both phi-segments.
 
          XS1  = R * COSPH1
          YS1  = R * SINPH1
          XS2  = R * COSPH2
          YS2  = R * SINPH2
 
C......  Get distances (squared) from given point to each endpoint.
 
          DISTS1 = (X(1) - XS1)**2  +  (X(2) - YS1)**2
          DISTS2 = (X(1) - XS2)**2  +  (X(2) - YS2)**2
 
C......  Get distance to that phi-segment whose endpoint
C......  is closest to the given point.
 
          IF (DISTS1 .LE. DISTS2) THEN
            SAFSEG = ABS(SINPH1 * X(1) - COSPH1 * X(2))
          ELSE
            SAFSEG = ABS(SINPH2 * X(1) - COSPH2 * X(2))
          ENDIF
 
  100   CONTINUE
 
 
        IF (SAFER .LE. 0.0)  THEN
 
C---------------------------------------------------------------------------
C
C......  Here we handle the case in which  SAFER < 0, i.e., the point is
C......  inside the polygon but outside the inscribed polycone. We must
C......  do an accurate calculation of "SAFER".
C
C---------------------------------------------------------------------------
 
          FACT = SAFZ2 * ZLENI
          RAD1 = PT2 + FACT * (PT4 - PT2)
          RAD2 = PT3 + FACT * (PT5 - PT3)
          RR1 = RAD1 * RAD1
          RR2 = RAD2 * RAD2
 
          IF(X(1).EQ.0.)THEN
             PHI   = ATAN2(X(2), 1.E-8)
          ELSE
             PHI   = ATAN2(X(2), X(1))
          ENDIF
          IF (PHI .LT. PHI1)  PHI = PHI + TWOPI
 
          DIST=0.
          IF (PHI .GE. PHI1  .AND.  PHI .LE. PHI2)  THEN
            PHIREL = PHI - PHI1
            NSECTR = INT(PHIREL / DELPHI) + 1
            PHICTR  = PHI1 + (2.0*NSECTR - 1.0) * DPHI2
            COSPHC = COS (PHICTR)
            SINPHC = SIN (PHICTR)
 
 
            IF (R2 .GE. RR2) THEN
              SR2 = RAD2
              DIST = ABS (COSPHC * X(1) + SINPHC * X(2) - SR2)
              FACTC  = (PT5 - PT3) * ZLENI
            ELSEIF (R2 .LE. RR1) THEN
              RIN  = RAD1 * CSDPH2
              SRIN = RIN
              DIST = ABS (COSPHC * X(1) + SINPHC * X(2) - SRIN)
              FACTC = (PT4 - PT2) * ZLENI
            ELSE
              PRINT *, 'GNPGON: FACTC not calculated!!!'
              FACTC=0
            ENDIF
 
          ENDIF
 
          SAFER  = DIST / SQRT(1.0 + FACTC*FACTC)
 
        ENDIF
 
        SAFE  = MIN (SAFEZ, SAFER, SAFSEG)
 
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (Nierhaus)  |
C     ------------------------------------------------
 
 
      CALL GNPGO1(X,PAR,SNXT)
C
  999 CONTINUE
      END
 
