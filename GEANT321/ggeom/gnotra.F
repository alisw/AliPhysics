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
      SUBROUTINE GNOTRA (X, PAR, IACT, IFLAG, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH 'TRD1' OR 'TRD2' *
C.    *       VOLUME, FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6) *
C.    *                                                                *
C.    *       PAR   (input)  : volume parameters                       *
C.    *       IACT  (input)  : action flag                             *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, and SNXT only if SNEXT .GT.new SAFE *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *       IFLAG (input)  : = 1 for TRD1, 2 for TRD2                *
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
      DIMENSION X(6),PAR(5)
C.
C.    ----------------------------------------------------------------
C.
 
      SNXT = BIG
      IF (IFLAG .EQ. 1) THEN
        FACTX  = (PAR(2) - PAR(1)) / (2.0*PAR(4))
      ELSEIF (IFLAG .EQ. 2) THEN
        FACTX  = (PAR(2) - PAR(1)) / (2.0*PAR(5))
        FACTY  = (PAR(4) - PAR(3)) / (2.0*PAR(5))
      ENDIF
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
        IF (IFLAG .EQ. 1) THEN
C*******************************************************
C
C......  Shape "TRD1":   only x varies with z.
C
C*******************************************************
 
          SAF2  = MAX (ABS(X(2)) - PAR(3),  0.0)
          SAF3  = MAX (ABS(X(3)) - PAR(4),  0.0)
C
C......  Distance along x-direction to slanted wall.
C
          DISTX = PAR(1) + FACTX * (PAR(4) + X(3))
          IF (ABS(X(1)) .GT. DISTX) THEN
            SAF1  = (ABS(X(1)) - DISTX) / SQRT(1.0 + FACTX*FACTX)
          ELSE
            SAF1  = 0.0
          ENDIF
 
 
        ELSE
C*******************************************************
C
C......  Shape "TRD2":  both x and y vary with z.
C
C*******************************************************
 
          SAF3  = MAX (ABS(X(3)) - PAR(5),  0.0)
 
C
C......  Distance along x-direction to slanted wall.
C
          DISTX = PAR(1) + FACTX * (PAR(5) + X(3))
          IF (ABS(X(1)) .GT. DISTX) THEN
            SAF1  = (ABS(X(1)) - DISTX) / SQRT(1.0 + FACTX*FACTX)
          ELSE
            SAF1  = 0.0
          ENDIF
C
C......  Distance along y-direction to slanted wall.
C
          DISTY = PAR(3) + FACTY * (PAR(5) + X(3))
          IF (ABS(X(2)) .GT. DISTY) THEN
            SAF2  = (ABS(X(2)) - DISTY) / SQRT(1.0 + FACTY*FACTY)
          ELSE
            SAF2  = 0.0
          ENDIF
 
        ENDIF
 
        SAFE  = MAX (SAF1, SAF2, SAF3)
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
C
C                 Only x varies with z.
C
      IF(IFLAG.NE.1) GO TO 50
      S3=PAR(4)
      FACTX = (PAR(2) - PAR(1)) / (2.0*S3)
      S1  = 0.5 * (PAR(1) + PAR(2))  +  FACTX * X(3)
      DS11=X(4) - FACTX * X(6)
      DS12=X(4) + FACTX * X(6)
      S2=PAR(3)
      DS21=X(5)
      DS22=X(5)
      GO TO 60
C
C            Both x and y widths vary with z.
C
   50 CONTINUE
      S3=PAR(5)
      FACTX = (PAR(2) - PAR(1)) / (2.0*S3)
      FACTY = (PAR(4) - PAR(3)) / (2.0*S3)
      S1 = 0.5 * (PAR(1) + PAR(2))  +  FACTX * X(3)
      DS11=X(4) - FACTX * X(6)
      DS12=X(4) + FACTX * X(6)
      S2 = 0.5 * (PAR(3) + PAR(4))  +  FACTY * X(3)
      DS21=X(5) - FACTY * X(6)
      DS22=X(5) + FACTY * X(6)
 
C
C            FIRST DO THE Z COORDINATE.
C
   60 SMAX=9.0E+10
      SMIN=0.0
      IF(X(6).EQ.0.0.AND.S3.LT.ABS(X(3))) GO TO 999
      IF(X(6).EQ.0.0) GO TO 70
      SMIN=(S3-X(3))/X(6)
      SMAX=(-S3-X(3))/X(6)
      IF(SMAX.GT.SMIN) GO TO 70
      SN=SMAX
      SMAX=SMIN
      SMIN=SN
C
   70 IF(SMIN.LT.0.0) SMIN=0.0
      IF(SMAX.LT.SMIN) GO TO 999
C
C             NOW DO X COORDINATE.
C
      SS1=S1-X(1)
      SS2=-S1-X(1)
      SN1=-1.0
      SN2=SMAX
*
      IF(SS1.LT.0.0) SN1=SMAX
      IF(SS2.GT.0.0) SN2=-1.0
*
      IF(DS11.NE.0.0) SN1=SS1/DS11
      IF(DS12.NE.0.0) SN2=SS2/DS12
      IF(DS12*DS11.EQ.0.0) GO TO 90
      IF(DS12.LT.0.0.AND.DS11.LT.0.0) GO TO 110
      ST=SN2
      SN2=SN1
      SN1=ST
      IF(DS12.GT.0.0.AND.DS11.GT.0.0) GO TO 110
      IF(DS12.LT.0.0) GO TO 80
      IF(SN2.GT.SN1) SN1=SN2
      SN2=SMAX
      GO TO 110
C
   80 CONTINUE
C
      IF(SN2.GT.SN1) SN2=SN1
      SN1=0.0
      GO TO 110
C
   90 CONTINUE
C
      IF(DS12.EQ.0.0.AND.SS2.GT.0.0) GO TO 999
      IF(DS12.EQ.0.0.AND.DS11.LT.0.0) GO TO 110
      IF(DS11.EQ.0.0.AND.SS1.LT.0.0) GO TO 999
      IF(DS11.EQ.0.0.AND.DS12.LE.0.0) GO TO 110
      IF(DS11.EQ.0.0) GO TO 100
      SN2=SN1
      SN1=0.0
      GO TO 110
C
  100 CONTINUE
      SN1=SN2
      SN2=SMAX
C
  110 CONTINUE
C
      IF(SN1.GT.SMIN) SMIN=SN1
      IF(SN2.LT.SMAX) SMAX=SN2
      IF(SMAX.LE.SMIN) GO TO 999
C
C             NOW DO Y COORDINATE.
C
      SS1=S2-X(2)
      SS2=-S2-X(2)
      SN1=-1.0
      SN2=SMAX
*
      IF(SS1.LT.0.0) SN1=SMAX
      IF(SS2.GT.0.0) SN2=-1.0
*
      IF(DS21.NE.0.0) SN1=SS1/DS21
      IF(DS22.NE.0.0) SN2=SS2/DS22
      IF(DS22*DS21.EQ.0.0) GO TO 130
      IF(DS22.LT.0.0.AND.DS21.LT.0.0) GO TO 150
      ST=SN2
      SN2=SN1
      SN1=ST
      IF(DS22.GT.0.0.AND.DS21.GT.0.0) GO TO 150
      IF(DS22.LT.0.0) GO TO 120
      IF(SN2.GT.SN1) SN1=SN2
      SN2=SMAX
      GO TO 150
C
  120 CONTINUE
C
      IF(SN2.GT.SN1) SN2=SN1
      SN1=0.0
      GO TO 150
C
  130 CONTINUE
C
      IF(DS22.EQ.0.0.AND.SS2.GT.0.0) GO TO 999
      IF(DS22.EQ.0.0.AND.DS21.LT.0.0) GO TO 150
      IF(DS21.EQ.0.0.AND.SS1.LT.0.0) GO TO 999
      IF(DS21.EQ.0.0.AND.DS22.LE.0.0) GO TO 150
      IF(DS21.EQ.0.0) GO TO 140
      SN2=SN1
      SN1=0.0
      GO TO 150
C
  140 CONTINUE
      SN1=SN2
      SN2=SMAX
C
  150 CONTINUE
C
      IF(SN1.GT.SMIN) SMIN=SN1
      IF(SN2.LT.SMAX) SMAX=SN2
C
      IF(SMAX.LE.SMIN) GO TO 999
      IF(SMIN.LE.0.)GO TO 999
      SNXT = SMIN
 
 999  CONTINUE
      END
