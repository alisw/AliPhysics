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
      SUBROUTINE GNOCTU(X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH 'CTUB'           *
C.    *      VOLUME, FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6)  *
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
C.    *    ==>Called by : GNEXT, GNOPCO, GTNEXT                        *
C.    *         Author   A.McPherson   ********                        *
C.    *       MODIFICATION LOG.                                        *
C.    *       18-July-89 M.Guckes modifications due to GEANG 3.13      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6),PAR(11),SNP(2),CSP(2)
*
*     ------------------------------------------------------------------
*
      SNXT = BIG
*
* *** Make sure that (-pi).LT.PHI.GE.(+pi)
*
      PH1  = MOD(PAR(4)*DEGRAD,TWOPI)
      PH2  = MOD(PAR(5)*DEGRAD,TWOPI)
      DPHI = PH2 - PH1
      IF(DPHI.LT.TWOPI) THEN
         DPHI = PH2 - PH1
         SNP(1) = SIN(PH1)
         CSP(1) = COS(PH1)
         SNP(2) = SIN(PH2)
         CSP(2) = COS(PH2)
      ELSE
         PH1  = 0.
         PH2  = TWOPI
         DPHI = TWOPI
         SNP(1) = 0.
         CSP(1) = 1.
         SNP(2) = 0.
         CSP(2) = 1.
      ENDIF
C
*
* *** Check that current point is outside the CTUB.  compute SAFE
*
      R02=X(1)*X(1)+X(2)*X(2)
      R=SQRT(R02)
      SAF1=PAR(1)-R
      SAF2=R-PAR(2)
      SAF3 = ( X(3)+PAR(3) )*PAR(8)+X(1)*PAR(6)+X(2)*PAR(7)
      SAF4 = ( X(3)-PAR(3) )*PAR(11)+X(1)*PAR(9)+X(2)*PAR(10)
*
      IF(IACT.EQ.3) GO TO 20
*
*     Phi segment (P. Weidhaas)
*       We compute here the distance (SAF5)
*        to the PHI-segment boundary that is closest to the point:
*
*        SAF5 is only calculated if PHI lies outside the interval
*        [PH1, PH2]. Here PHI is the angle to the given point
*        (thus we only consider SAF5 if the point is outside the
*        PHI-segment : FIOUT > 0.
*        Algorithm to find SAF5 (same as in routine "GNTUBE"):
*        For each PHI-boundary we find the distance from the given
*        point to the outer (at RMAX) point of the segment boundary
*        (DSP1 and DSP2, resp.). If DSP1 < DSP2, we define
*        SAF5 to be the distance to segment PHI1; else we set
*        SAF5 to be the distance to segment PHI2.
      SAF5=0.
      IF(R.GT.0.) THEN
         PHI    = ATAN2(X(2),X(1))
         FIOUT  = DPHI*(PHI-PH1)*(PHI-PH2)
         IF(FIOUT.LE.0.) GO TO 16
         DSP1 = (X(1)-PAR(2)*CSP(1))**2 + (X(2)-PAR(2)*SNP(1))**2
         DSP2 = (X(1)-PAR(2)*CSP(2))**2 + (X(2)-PAR(2)*SNP(2))**2
         IF(DSP1.LE.DSP2) THEN
            SAF5 = ABS(X(1)*SNP(1) - X(2)*CSP(1))
         ELSE
            SAF5 = ABS(X(1)*SNP(2) - X(2)*CSP(2))
         ENDIF
      ENDIF
   16 CONTINUE
*
*     Compute SAFE
      SAFE = 0.
      IF(SAF1.GT.SAFE) SAFE = SAF1
      IF(SAF2.GT.SAFE) SAFE = SAF2
      IF(SAF3.GT.SAFE) SAFE = SAF3
      IF(SAF4.GT.SAFE) SAFE = SAF4
      IF(SAF5.GT.SAFE) SAFE = SAF5
*
*     Point inside the volume ?
      IF(SAFE.EQ.0.) THEN
         SNXT = -10.
         GO TO 999
      ENDIF
*
      IF(IACT.EQ.0)                   GO TO 999
      IF(IACT.EQ.1.AND.SAFE.GT.SNEXT) GO TO 999
*
   20 CONTINUE
*
* *** Compute SNXT
*
      SMAX=BIG
      V1 = X(4)*PAR(6)+X(5)*PAR(7)+X(6)*PAR(8)
      V2 = X(4)*PAR(9)+X(5)*PAR(10)+X(6)*PAR(11)
      IF( V1 .GE. 0.0 .AND. SAF3 .GE. 0.0 ) GO TO 999
      IF( V2 .GE. 0.0 .AND. SAF4 .GE. 0.0 ) GO TO 999
C
      IF( SAF3 .GT. 0.0 .AND. SAF4 .LE. 0.0 ) THEN
        SMIN = -SAF3/V1
        IF( V2 .GT. 0 ) THEN
          SMAX = -SAF4/V2
        ELSE
          SMAX = BIG
        ENDIF
      ELSEIF( SAF4 .GT. 0.0 .AND. SAF3 .LE. 0.0 ) THEN
        SMIN = -SAF4/V2
        IF( V1 .GT. 0 ) THEN
          SMAX = -SAF3/V1
        ELSE
          SMAX = BIG
        ENDIF
      ELSEIF( SAF3 .LE. 0.0 .AND. SAF4 .LE. 0.0 ) THEN
        SMIN = 0.0
        SMAX = BIG
        IF( V1 .GT. 0.0 ) SMAX = -SAF3/V1
        IF( V2 .GT. 0.0 .AND. SMAX .GT. -SAF4/V2 ) SMAX = -SAF4/V2
      ENDIF
*
      SMIN1=SMIN
      SMAX1=-1.0
      SMIN2=SMIN
      SMAX2=SMAX
C
      DXY2=(1+X(6))*(1-X(6))
      IF(DXY2.LT.1.0E-10.AND.(R.LT.PAR(1)
     +   .OR.R.GT.PAR(2)))GO TO 999
      IF(DXY2.EQ.0.) GO TO 30
C
      BA=(X(1)*X(4)+X(2)*X(5))/DXY2
      IF(R.GE.PAR(2).AND.BA.GE.0.0) GO TO 999
C
      CA=(R02-PAR(2)*PAR(2))/DXY2
      DISC=BA*BA-CA
      IF(DISC.LE.0.0) GO TO 999
C
      RDISC=SQRT(DISC)
      SMAR=-BA+RDISC
      IF(SMAR.LT.SMAX) SMAX=SMAR
      SMIR=-BA-RDISC
      IF(SMIR.GT.SMIN) SMIN=SMIR
      IF(SMAX.LE.SMIN) GO TO 999
C
      CA=(R02-PAR(1)*PAR(1))/DXY2
      DISC=BA*BA-CA
C
      SMIN1=SMIN
      SMAX1=-1.0
      SMIN2=SMIN
      SMAX2=SMAX
C
      IF(DISC.LE.0.0) GO TO 30
      RDISC=SQRT(DISC)
      SMI2=-BA+RDISC
      IF(SMI2.GT.SMIN2)SMIN2=SMI2
      SMAX1=-BA-RDISC
      IF(SMAX.LT.SMAX1) SMAX1=SMAX
C
   30 CONTINUE
C
      IF( PAR(4) .EQ. 0.0 .AND. PAR(5) .EQ. 360.0 ) GO TO 120
C
C                   Now do the phi stuff.
C
      DPSGN=X(1)*X(5)-X(2)*X(4)
      ISMIN=0
      ISMAX=0
C
      TSGN=SNP(1)
      TCSG=CSP(1)
      DEN=X(4)*TSGN-X(5)*TCSG
      IF(DEN.EQ.0.0) GO TO 40
      SNL=(X(2)*TCSG-X(1)*TSGN)/DEN
      IF(ABS(TSGN).GT.1.E-6.AND.(X(2)+SNL*X(5))*TSGN.LT.0.) GO TO 40
      IF(ABS(TCSG).GT.1.E-6.AND.(X(1)+SNL*X(4))*TCSG.LT.0.) GO TO 40
C
      ISMIN=1
      SMIN=SNL
      IF(DPSGN.GT.0.0) GO TO 40
      ISMIN=0
      ISMAX=1
      SMAX=SNL
C
   40 CONTINUE
C
      TSGN=SNP(2)
      TCSG=CSP(2)
      DEN=X(4)*TSGN-X(5)*TCSG
      IF(DEN.EQ.0.0) GO TO 60
      SNH=(X(2)*TCSG-X(1)*TSGN)/DEN
      IF(ABS(TSGN).GT.1.E-6.AND.(X(2)+SNH*X(5))*TSGN.LT.0.) GO TO 60
      IF(ABS(TCSG).GT.1.E-6.AND.(X(1)+SNH*X(4))*TCSG.LT.0.) GO TO 60
      IF(DPSGN.LT.0.0) GO TO 50
      ISMAX=1
      SMAX=SNH
      GO TO 60
C
   50 CONTINUE
      ISMIN=1
      SMIN=SNH
C
   60 CONTINUE
C
      IF(ISMIN.EQ.0.OR.ISMAX.EQ.0) GO TO 90
      IF(SMAX.LT.0.0.AND.SMAX.GT.SMIN) GO TO 999
      IF(SMIN.LT.0.0) SMIN=0.0
      IF(SMAX.LT.0.0) GO TO 110
      IF(SMAX.GT.SMIN) GO TO 70
C
C                         SMAX +VE and less than SMIN - 2 allowed
C                         regions.
C
      IF(SMAX.GT.SMIN1.AND.SMAX1.GT.SMIN1) GO TO 130
      IF(SMIN.GT.SMIN1) SMIN1=SMIN
      IF(SMAX1.GT.SMIN1) GO TO 130
      IF(SMAX.GT.SMIN2) GO TO 140
      IF(SMIN.GT.SMAX2) GO TO 999
      IF(SMIN.GT.SMIN2) SMIN2=SMIN
      GO TO 140
C
   70 CONTINUE
      IF(SMIN.GT.SMIN1) SMIN1=SMIN
      IF(SMAX.LT.SMAX1) SMAX1=SMAX
      IF(SMAX1.LT.SMIN1) GO TO 80
      GO TO 130
C
   80 CONTINUE
      IF(SMIN.GT.SMIN2) SMIN2=SMIN
      IF(SMAX.LT.SMAX2) SMAX2=SMAX
      IF(SMAX2.LT.SMIN2) GO TO 999
      GO TO 140
C
   90 CONTINUE
      IF(ISMIN.EQ.1) GO TO 110
      IF(ISMAX.EQ.0) GO TO 100
      IF(SMAX.LT.SMAX1) SMAX1=SMAX
      IF(SMAX.LT.SMAX2) SMAX2=SMAX
      GO TO 120
C
  100 CONTINUE
      DPH=PAR(5)-PAR(4)
      IF(DPH.LT.180.0.AND.DPH.GT.0.0) GO TO 999
      IF(DPH.LT.-180.0) GO TO 999
      GO TO 120
C
  110 CONTINUE
      IF(SMIN.GT.SMIN1) SMIN1=SMIN
      IF(SMAX1.GT.SMIN1) GO TO 130
      IF(SMIN.GT.SMIN2) SMIN2=SMIN
      IF(SMAX2.LT.SMIN2) GO TO 999
      GO TO 140
C
  120 CONTINUE
      IF(SMAX1.GT.SMIN1) GO TO 130
      IF(SMAX2.GT.SMIN2) GO TO 140
      GO TO 999
C
  130 CONTINUE
      IF(SMIN1.LE.0.)GO TO 999
      IF(SNXT.GT.SMIN1) SNXT=SMIN1
      GO TO 999
C
  140 CONTINUE
      IF(SMIN2.LE.0.)GO TO 999
      IF(SNXT.GT.SMIN2) SNXT=SMIN2
C
  999 CONTINUE
      END
