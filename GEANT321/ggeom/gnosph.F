*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNOSPH (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'SPHE' VOLUME,  *
C.    *       FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6)SPHERE   *
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
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#include "geant321/gconsp.inc"
      REAL X(6),PAR(6),SNEXT,SNXT,SAFE
      DIMENSION IS(4),SS(4)
C.
C.    ----------------------------------------------------------------
C.
 
      SNXT = BIG
      R2  = X(1)*X(1) + X(2)*X(2) + X(3)*X(3)
      R   = SQRT (R2)
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
        RIN  = PAR(1)
        ROUT = PAR(2)
        IF (R .LT. RIN) THEN
          SAFE  = RIN - R
        ELSEIF (R .GT. ROUT) THEN
          SAFE  = R - ROUT
        ENDIF
 
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
      BA=X(1)*X(4)+X(2)*X(5)+X(3)*X(6)
      IF(R2.GE.PAR(2)*PAR(2).AND.BA.GE.0.0) GO TO 999
C
      CA=R2-PAR(2)*PAR(2)
      DISC=BA*BA-CA
      IF(DISC.LE.0.0) GO TO 999
C
      RDISC=SQRT(DISC)
      SMAX=-BA+RDISC
      SMIN=-BA-RDISC
C
C                   NOW DO RMIN
C
      CA=R2-PAR(1)*PAR(1)
      DISC=BA*BA-CA
C
      SMIN1=SMIN
      SMAX1=-1.0
      SMIN2=SMIN
      SMAX2=SMAX
C
      IF(DISC.LE.0.0) GO TO 30
      RDISC=SQRT(DISC)
      SMIN2=-BA+RDISC
      SMAX1=-BA-RDISC
C
   30 CONTINUE
C
C                   NOW DO THE PHI STUFF.
C
      IP2=0
      SMNP1=0.0
      SMXP1=SMAX2
C
      IF(PAR(6)-PAR(5).GE.360.0) GO TO  110
C
      DPSGN=X(1)*X(5)-X(2)*X(4)
      PHL=PAR(5)/RADDEG
      PHH=PAR(6)/RADDEG
      ISMIN=0
      ISMAX=0
C
      TSGN=SIN(PHL)
      TCSG=COS(PHL)
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
      TSGN=SIN(PHH)
      TCSG=COS(PHH)
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
      IF(ISMIN.EQ.0.OR.ISMAX.EQ.0) GO TO 80
      IF(SMAX.LT.0.0.AND.SMAX.GT.SMIN) GO TO 999
      IF(SMIN.LT.0.0) SMIN=0.0
      IF(SMAX.LT.0.0) GO TO  100
      IF(SMAX.GT.SMIN) GO TO 70
C
C                         SMAX +VE AND LESS THAN SMIN - 2 ALLOWED
C                         REGIONS.
C
      IP2=1
      SMXP1=SMAX
      SMNP2=SMIN
      SMXP2=SMAX2
      GO TO  110
C
   70 CONTINUE
C
C                         SMIN +VE AND SMAX GT SMIN: NORMAL SINGLE
C                         REGION
C
      SMNP1=SMIN
      SMXP1=SMAX
      GO TO  110
C
   80 CONTINUE
      IF(ISMIN.EQ.1) GO TO  100
      IF(ISMAX.EQ.0) GO TO  90
C
C                    SMAX BUT NO SMIN
C
      SMXP1=SMAX
      GO TO  110
C
   90 CONTINUE
C
C                   NO SMIN OR SMAX: ALWAYS IN OR ALWAYS OUT.
C
      DPH=PAR(5)-PAR(4)
      IF(DPH.LT.180.0.AND.DPH.GT.0.0) GO TO 999
      IF(DPH.LT.-180.0) GO TO 999
      GO TO  110
C
  100 CONTINUE
C
C                  SMIN BUT NO SMAX.
C
      SMNP1=SMIN
C
  110 CONTINUE
C
C                  NOW DO THETA.
C
      IT2=0
      SMNT1=0.0
      SMXT1=SMAX2
      IF(PAR(4)-PAR(3).GE.180.0) GO TO 360
C
      TH=PAR(3)
      IT=1
      ITLN=0
      ITLX=0
      ITHN=0
      ITHX=0
C
  120 CONTINUE
C
      IF(TH.NE.90.0) GO TO 130
      IF(X(6).EQ.0.0) GO TO 220
C
      ST=-X(3)/X(6)
      STST=-X(6)
      GO TO 180
C
  130 CONTINUE
C
      TT=TAN(TH/RADDEG)
      TT2=TT*TT
C
      A=X(4)*X(4)+X(5)*X(5)-TT2*X(6)*X(6)
      B=X(1)*X(4)+X(2)*X(5)-TT2*X(3)*X(6)
      C=X(1)*X(1)+X(2)*X(2)-TT2*X(3)*X(3)
C
      IF(A.NE.0.0) GO TO 140
      IF(B.EQ.0.0) GO TO 220
C
      ST=-C*0.5/B
C
      Z=X(3)+ST*X(6)
      IF(Z*TT.LT.0.0) GO TO 220
C
      STST=(B+ST*A)/Z
      ITRY=2
C
      GO TO 180
  140 CONTINUE
C
      BA=B/A
      CA=C/A
      DISC=BA*BA-CA
      IF(DISC.LT.0.0) GO TO 220
C
      RDISC=0.0
      IF(DISC.GT.0.0) RDISC=SQRT(DISC)
      ITRY=1
      ST=-BA-RDISC
C
  150 CONTINUE
C
      IF(ST.LT.0.0) GO TO 160
      Z=X(3)+ST*X(6)
      IF(Z.EQ.0.0.AND.ABS(A).LT.0.0) GO TO 170
      IF(RDISC.EQ.0.0) GO TO 160
      IF(Z*TT.LT.0.0) GO TO 160
C
      STST=(B+ST*A)/Z
      GO TO 180
C
  160 CONTINUE
C
      IF(ITRY.EQ.2) GO TO 220
      ST=-BA+RDISC
      ITRY=2
      GO TO 150
C
  170 CONTINUE
      STST=-X(6)
  180 CONTINUE
C
      IF(IT.NE.1) GO TO 200
      IF(STST.GT.0.0) GO TO 190
      ITLX=1
      SMXTL=ST
      GO TO 160
C
  190 CONTINUE
      ITLN=1
      SMNTL=ST
      GO TO 160
C
  200 CONTINUE
      IF(STST.GT.0.0) GO TO 210
      ITHN=1
      SMNTH=ST
      GO TO 160
C
  210 CONTINUE
      ITHX=1
      SMXTH=ST
      GO TO 160
C
  220 CONTINUE
      IF(IT.EQ.2) GO TO 230
      IT=2
      TH=PAR(4)
      GO TO 120
C
  230 CONTINUE
C
C              ORDER THE VARIOUS BOUNDARIES.
C
      ICOUNT=0
      IF(ITLN.EQ.0.OR.SMNTL.LE.0.0) GO TO 240
      IS(1)=1
      SS(1)=SMNTL
      ICOUNT=1
C
  240 CONTINUE
      IF(ITLX.EQ.0.OR.SMXTL.LE.0.0) GO TO 260
      IPL=ICOUNT+1
      IF(ICOUNT.EQ.0.OR.SMXTL.GT.SS(1)) GO TO 250
      IS(2)=IS(1)
      SS(2)=SS(1)
      IPL=1
  250 CONTINUE
      ICOUNT=ICOUNT+1
      IS(IPL)=2
      SS(IPL)=SMXTL
C
  260 CONTINUE
      IST=3
      IF(ITHN.EQ.0.OR.SMNTH.LE.0.0) GO TO 320
      STEST=SMNTH
C
  270 CONTINUE
      IPL=ICOUNT+1
      IF(ICOUNT.EQ.0) GO TO 310
      DO 280 IC=1,ICOUNT
      IC1=ICOUNT-IC+1
      IF(STEST.GT.SS(IC1)) GO TO 290
      IPL=IPL-1
  280 CONTINUE
C
  290 CONTINUE
      IF(IPL.EQ.ICOUNT+1) GO TO 310
      IM=ICOUNT+1-IPL
      DO 300 I=1,IM
      I1=ICOUNT-I+1
      I2=I1+1
      IS(I2)=IS(I1)
      SS(I2)=SS(I1)
  300 CONTINUE
C
  310 CONTINUE
      ICOUNT=ICOUNT+1
      IS(IPL)=IST
      SS(IPL)=STEST
C
  320 CONTINUE
      IF(IST.EQ.4) GO TO 330
      IF(ITHX.EQ.0.OR.SMXTH.LE.0.0) GO TO 330
      IST=4
      STEST=SMXTH
      GO TO 270
C
  330 CONTINUE
C
C               CHECK WHETHER 1ST IS MAX OR MIN.
C
      IF(ICOUNT.EQ.0) GO TO 350
      IF(IS(1).EQ.2.OR.IS(1).EQ.4) GO TO 340
C
C               START WITH MIN.
C
      SMNT1=SS(1)
      IF(ICOUNT.GE.2) SMXT1=SS(2)
      IF(ICOUNT.LE.2) GO TO 360
      IT2=1
      SMNT2=SS(3)
      SMXT2=SMAX2
      IF(ICOUNT.GE.4) SMXT2=SS(4)
      GO TO 360
C
  340 CONTINUE
C
C              START WITH MAX SO 1ST MIN IS 0.0
C
      SMNT1=0.0
      SMXT1=SS(1)
      IF(ICOUNT.LE.1) GO TO 360
      IT2=1
      SMNT2=SS(2)
      SMXT2=SMAX2
      IF(ICOUNT.GE.3) SMXT2=SS(3)
      GO TO 360
C
  350 CONTINUE
C
C              NO INTERSECTIONS ALWAYS IN OR ALWAYS OUT.
C
      R=X(1)*X(1)+X(2)*X(2)
      IF(R.GT.0.0) R=SQRT(R)
      TH=90.0
      IF(X(3).NE.0.0) TH=ATAN(R/X(3))*RADDEG
      IF(TH.LT.0.0) TH=180.0+TH
      IF(TH.LT.PAR(3).OR.TH.GT.PAR(4)) GO TO 999
C
  360 CONTINUE
C
C              NOW FIND SMALLEST S ALOWED BY ALL.
C
      IF(SMAX1.LE.SMIN1) GO TO 370
      SMAXR=SMAX1
      SMINR=SMIN1
      IRT=1
      GO TO 380
C
  370 CONTINUE
      SMAXR=SMAX2
      SMINR=SMIN2
      IRT=2
C
  380 CONTINUE
      IF(SMNP1.GT.SMAXR) GO TO 430
      IF(SMXP1.LT.SMINR) GO TO 390
      SMIN=SMINR
      SMAX=SMAXR
      IF(SMNP1.GT.SMIN) SMIN=SMNP1
      IF(SMXP1.LT.SMAX) SMAX=SMXP1
      IPT=1
      GO TO 400
C
  390 CONTINUE
      IF(IP2.EQ.0) GO TO 430
      IF(SMNP2.GT.SMAXR) GO TO 430
      IF(SMXP2.LT.SMINR) GO TO 430
      SMIN=SMINR
      SMAX=SMAXR
      IF(SMNP2.GT.SMIN) SMIN=SMNP2
      IF(SMXP2.LT.SMAX) SMAX=SMXP2
      IPT=2
C
  400 CONTINUE
C
      IF(SMNT1.GT.SMAX) GO TO 420
      IF(SMXT1.LT.SMIN) GO TO 410
      IF(SMNT1.GT.SMIN) SMIN=SMNT1
      GO TO 440
C
  410 CONTINUE
      IF(IT2.EQ.0) GO TO 420
      IF(SMNT2.GT.SMAX) GO TO 420
      IF(SMXT2.LT.SMIN) GO TO 420
      IF(SMNT2.GT.SMIN) SMIN=SMNT2
      GO TO 440
C
  420 CONTINUE
      IF(IPT.EQ.1) GO TO 390
  430 CONTINUE
      IF(IRT.EQ.1) GO TO 370
      GO TO 999
C
  440 CONTINUE
      IF(SMIN.LE.0.)GO TO 999
      SNXT = SMIN
 
  999 CONTINUE
      END
