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
      SUBROUTINE GNOPCO (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'PCON' VOLUME,  *
C.    *       FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6)         *
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
      DIMENSION X(6), PAR(9), PT(7), XT(6)
 
      EQUIVALENCE (PT(1),PT1), (PT(2),PT2), (PT(3),PT3), (PT(4),PT4)
      EQUIVALENCE (PT(5),PT5), (PT(6),PT6), (PT(7),PT7)
C.
C.   --------------------------------------------------------
C.
 
      SNXT = BIG
      SAFE = BIG
      R2 = X(1)*X(1) + X(2)*X(2)
      R  = SQRT (R2)
      NZ  = PAR(3)
      NSEC= NZ - 1
      PT6=PAR(1)
      PT7=PAR(1)+PAR(2)
 
      IFLAG = 2
      IF(PAR(2).EQ.360.0) IFLAG = 1
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
 
        SAFEZ  = 0.0
C
C......  Obtain axial distance "SAFEZ".
C
        ZMIN   = PAR(4)
        ZMAX   = PAR(3*NZ+1)
 
        IF (X(3) .LT. ZMIN) THEN
          SAFEZ  = ZMIN - X(3)
        ELSEIF (X(3) .GT. ZMAX) THEN
          SAFEZ  = X(3) - ZMAX
        ENDIF
 
 
C----------------------------------------------------
C......  Prepare parameters for PHI-segmented cone:
C----------------------------------------------------
 
      IF (IFLAG .EQ. 2) THEN
 
        PHI1 = MOD (PT6+360.0 , 360.0)
        PHI2 = MOD (PT7+360.0 , 360.0)
        IF ( X(1).NE.0.0 .OR. X(2).NE.0.0 ) THEN
           PHI = ATAN2( X(2), X(1) ) * RADDEG
        ELSE
           PHI = 0.0
        ENDIF
        PHI  = MOD (PHI+360.0 , 360.0)
 
        SINPH1 = SIN(PHI1*DEGRAD)
        COSPH1 = COS(PHI1*DEGRAD)
        SINPH2 = SIN(PHI2*DEGRAD)
        COSPH2 = COS(PHI2*DEGRAD)
 
 
C......  Set flag "IOPENP" if point (X(1),X(2),X(3)) lies in the
C......  PHI-opening.
 
        IOPENP = 0
        IF (PHI2 .GT. PHI1) THEN
          IF (PHI.GT.PHI2 .OR. PHI.LT.PHI1) IOPENP = 1
        ELSE
          IF (PHI.GT.PHI2 .AND. PHI.LT.PHI1) IOPENP = 1
        ENDIF
      ENDIF
 
 
 
C------------------   Start of loop over Z-sections  --------------
 
        IPZ2=4
 
        DO 150 IS=1,NSEC
 
        IPZ1=IPZ2
        IPZ2=IPZ1+3
 
        SAFSEG = 0.0
        SAFER  = 0.0
 
        XT3=X(3) - 0.5 * (PAR(IPZ1)+PAR(IPZ2))
        DZ   = 0.5 * (PAR(IPZ2)-PAR(IPZ1))
 
        PT2  = PAR(IPZ1+1)
        PT3  = PAR(IPZ1+2)
        PT4  = PAR(IPZ2+1)
        PT5  = PAR(IPZ2+2)
C**** check DZ=0 segments
        IF (DZ.LE.0.) THEN
           IF ((R-PT2)*(R-PT4).LE.0. .OR. (R-PT3)*(R-PT5).LE.0.) THEN
              SAFER = ABS(XT3)
           ELSE
              GO TO 150
           ENDIF
           GO TO 100
        ENDIF
 
        IF (PT2 .NE. PT4) GO TO 50
        IF (PT3 .NE. PT5) GO TO 50
 
C**********************************************************
C
C......  The segment is a tube;  invoke the algorithm
C......  from "GNOTUB" inline to get "SAFER" and "SAFSEG".
C
C**********************************************************
 
        IF (R .LT. PT2) SAFER = PT2 - R
        IF (R .GT. PT3) SAFER = R - PT3
 
        IF (IFLAG .EQ. 2)  THEN
 
C********************************************************************
C......  Handle the case in which we have a PHI-segment of a tube.
C......  In addition to the radial distance (SAFER) and the
C......  axial distance (SAFEZ) we compute here the distance (SAFSEG)
C......  to the PHI-segment boundary that is closest to the point:
C
C......  SAFSEG is only calculated if PHI lies outside the interval
C......  [PHI1, PHI2]. Here PHI is the angle to the given point
C......  (thus we only consider SAFSEG if the point is outside the
C......  PHI-segment).
C
C......  Algorithm to find SAFSEG (same as in routine "GNTUBE"):
C
C......  For each PHI-boundary we find the distance from the given
C......  point to the outer (at RMAX) point of the segment boundary
C......  (DISTS1 and DISTS2, resp.). If DISTS1 < DISTS2, we define
C......  SAFSEG to be the distance to segment PHI1; else we set
C......  SAFSEG to be the distance to segment PHI2.
C********************************************************************
 
 
C......  Next eliminate those points whose angle PHI places them
C......  inside the given PHI-segment (IOPENP = 0).
 
          IF (IOPENP .EQ. 0) GO TO 100
 
 
C......  Get coordinates of outer endpoints (at RMAX) of both PHI-segments.
 
          XS1  = PT3 * COSPH1
          YS1  = PT3 * SINPH1
          XS2  = PT3 * COSPH2
          YS2  = PT3 * SINPH2
 
C......  Get distances (squared) from the given point to each endpoint.
 
          DISTS1 = (X(1) - XS1)**2  +  (X(2) - YS1)**2
          DISTS2 = (X(1) - XS2)**2  +  (X(2) - YS2)**2
 
C......  Get distance to that PHI-segment whose endpoint
C......  is closest to the given point.
 
          IF (DISTS1 .LE. DISTS2) THEN
            SAFSEG = ABS(X(1) * SINPH1 - X(2) * COSPH1)
          ELSE
            SAFSEG = ABS(X(1) * SINPH2 - X(2) * COSPH2)
          ENDIF
 
        ENDIF
 
        GO TO 100
 
 
   50   CONTINUE
 
C*********************************************************
C
C......  The segment is a cone; invoke the algorithm
C......  from "GNOCON" inline to get "SAFER" and "SAFSEG".
C
C*********************************************************
 
        ZLENI = 0.5 / DZ
 
        FACT1  = (PT4 - PT2) * ZLENI
        FACT2  = (PT5 - PT3) * ZLENI
        RIN   = PT2 + FACT1 * (DZ + XT3)
        ROUT  = PT3 + FACT2 * (DZ + XT3)
 
        IF (R .LT. RIN) THEN
          SAFER = (RIN - R) / SQRT(1.0 + FACT1 * FACT1)
        ELSE
          IF (R .GT. ROUT) THEN
            SAFER = (R - ROUT) / SQRT(1.0 + FACT2 * FACT2)
          ENDIF
        ENDIF
 
 
        IF (IFLAG .EQ. 2)  THEN
 
C********************************************************************
C......  Handle the case in which we have a PHI-segment of a cone.
C......  In addition to the radial distance (SAFER) and the
C......  axial distance (SAFEZ) we compute here the distance (SAFSEG)
C......  to the PHI-segment boundary that is closest to the point:
C
C......  SAFSEG is only calculated if PHI lies outside the interval
C......  [PHI1, PHI2]. Here PHI is the angle to the given point
C......  (thus we only consider SAFSEG if the point is outside the
C......  PHI-segment).
C
C......  Algorithm to find SAFSEG (same as in routine "GNTUBE"):
C
C......  For each PHI-boundary we find the distance from the given
C......  point to the outer (at ROUT) point of the segment boundary
C......  (DISTS1 and DISTS2, resp.). If DISTS1 < DISTS2, we define
C......  SAFSEG to be the distance to segment PHI1; else we set
C......  SAFSEG to be the distance to segment PHI2.
C********************************************************************
 
 
C......  Next eliminate those points whose angle PHI places them
C......  inside the given PHI-segment (IOPENP = 0).
 
           IF (IOPENP .EQ. 0) GO TO 100
 
C......  Get coordinates of outer endpoints (at ROUT) of both PHI-segments.
 
          IF (XT3 .LT. -DZ) THEN
            ROUT = PT3
          ELSEIF (XT3 .GT. DZ) THEN
            ROUT = PT5
          ENDIF
 
          XS1  = ROUT * COSPH1
          YS1  = ROUT * SINPH1
          XS2  = ROUT * COSPH2
          YS2  = ROUT * SINPH2
 
C......  Get distances (squared) from the given point to each endpoint.
 
          DISTS1 = (X(1) - XS1)**2  +  (X(2) - YS1)**2
          DISTS2 = (X(1) - XS2)**2  +  (X(2) - YS2)**2
 
C......  Obtain distance to that PHI-segment whose endpoint
C......  is closest to the given point.
 
          IF (DISTS1 .LE. DISTS2) THEN
            SAFSEG = ABS(X(1) * SINPH1 - X(2) * COSPH1)
          ELSE
            SAFSEG = ABS(X(1) * SINPH2 - X(2) * COSPH2)
          ENDIF
 
        ENDIF
 
 
  100   CONTINUE
        TSAFE  = MAX (SAFEZ, SAFER, SAFSEG)
 
        IF (TSAFE .GT. 0.0) THEN
          IF (TSAFE .LT. SAFE) SAFE  = TSAFE
        ENDIF
 
        IF (TSAFE .EQ. 0.0) THEN
          IF (X(3) .LT. PAR(IPZ1)) GO TO 200
        ENDIF
  150   CONTINUE
 
  200   CONTINUE
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
      TSNXT = BIG
 
      DO 210 I=1, 6
        XT(I) = X(I)
  210 CONTINUE
 
      IPZ2 = 4
 
      DO 300  IS=1, NSEC
        IPZ1  = IPZ2
        IPZ2  = IPZ1 + 3
 
        XT(3) = X(3) - 0.5 * (PAR(IPZ1) + PAR(IPZ2))
 
        PT1 = 0.5 * (PAR(IPZ2) - PAR(IPZ1))
        IF (PT1 .LE. 0.0) GO TO 300
 
        IF (PAR(IPZ1+1) .NE. PAR(IPZ2+1)) GO TO 250
        IF (PAR(IPZ1+2) .NE. PAR(IPZ2+2)) GO TO 250
 
C
C......  This Z-section is a tube.
C
        PT3 = PT1
        PT1 = PAR(IPZ1+1)
        PT2 = PAR(IPZ1+2)
        PT4 = PT6
        PT5 = PT7
 
        CALL GNOTUB (XT, PT, 3, IFLAG, SNEXT, TSNXT, TSAFE)
        GO TO 280
 
 
  250   CONTINUE
C
C......   This Z-section is a cone.
C
        PT2 = PAR(IPZ1+1)
        PT3 = PAR(IPZ1+2)
        PT4 = PAR(IPZ2+1)
        PT5 = PAR(IPZ2+2)
 
        CALL GNOCON (XT, PT, 3, IFLAG, SNEXT, TSNXT, TSAFE)
 
  280 CONTINUE
      IF (TSNXT .LT. SNXT)  SNXT = TSNXT
  300 CONTINUE
 
  999 CONTINUE
      END
