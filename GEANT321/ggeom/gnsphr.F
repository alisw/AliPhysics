*
* $Id$
*
* $Log$
* Revision 1.3  1996/03/28 08:50:23  cernlib
* In the call to function MAX  use NULL istead of 0. to get the proper
* type of argument.
*
* Revision 1.2  1996/02/27 10:12:06  ravndal
* Precision problem (neg. Sqrt) solved
*
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNSPHR (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'SPHE' VOLUME,  *
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
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#include "geant321/gconsp.inc"
      REAL X(6), PAR(6), SNXT, SNEXT, SAFE
      PARAMETER (ZERO=0,ONE=1,SMALL=ONE/BIG)
 
C---------------------------------------------------------------
 
      SNXT = BIG
      R2   = X(1)*X(1) + X(2)*X(2) + X(3)*X(3)
      R    = SQRT (R2)
      RIN  = PAR(1)
      ROUT = PAR(2)
      SAF1 = R - RIN
      SAF2 = ROUT - R
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
         SAFE = MIN (SAF1, SAF2)
         IF( RIN .EQ. 0.0 ) THEN
            SAFE = ROUT-R
            SAF3 = SAFE
            SAF4 = SAFE
         ELSE
            SAF3 = R
            SAF4 = R
         ENDIF
         IF( R .GT. 1.E-5) THEN
            IF( PAR(4)-PAR(3).GE.180.) THEN
               SAF3 = BIG
            ELSE
               TH = ACOS(X(3)/SNGL(R))*RADDEG
               IF (TH.LT.0) TH = 180.+TH
               DTH1 = TH-PAR(3)
               DTH2 = PAR(4)-TH
               DTH = MIN(DTH1,DTH2)
               SAF3 = R*SIN(DTH*DEGRAD)
            ENDIF
            RXY2 = X(1)*X(1)+X(2)*X(2)
            IF( RXY2 .GT. 1.E-12.AND. PAR(6)-PAR(5) .LT. 360) THEN
               RXY = SQRT(RXY2)
               PHI=ATAN2(X(2),X(1))*RADDEG
               DPH = PHI-PAR(5)
               SG = SIGN(ONE,DPH)
               DPH = MOD( ABS(DPH), ONE*360 )
               IF(SG.LE.0.) DPH = 360.-DPH
               DPHT = PAR(6)-PAR(5)
               IF(DPHT .LT. 0.) DPHT = DPHT+360.
               DDPH = MIN(DPH,DPHT-DPH)
               IF( DDPH .GT. 90.) THEN
                  SAF4 = BIG
               ELSE
                  SAF4 = RXY*SIN(DDPH*DEGRAD)
               ENDIF
            ENDIF
         ENDIF
         SAFE = MIN(ONE*SAFE,SAF3,SAF4)
 
         IF (IACT .EQ. 0) GO TO 999
         IF (IACT .EQ. 1) THEN
            IF (SNEXT .LT. SAFE) GO TO 999
         ENDIF
      ENDIF
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
      IF(R.LT.1.0E-5) GO TO 70
 
      C1=(X(1)*X(4)+X(2)*X(5)+X(3)*X(6)) / R
      RTMP = PAR(2)
      SGN=ONE
      AC=R*R * (ONE-C1*C1)
      IF(AC.GT.PAR(1)**2.OR.C1.GT.0.0) GO TO 10
      SGN=-ONE
      RTMP = PAR(1)
   10 CONTINUE
      RTMPAC=MAX(RTMP*RTMP-AC,ZERO)
      SNXT=SQRT(RTMPAC)*SGN-R*C1
      DPSGN=X(1)*X(5)-X(2)*X(4)
      PHI2=PAR(6)
      IF(DPSGN.LT.0.0) PHI2=PAR(5)
C
      TSGN=SIN(PHI2*DEGRAD)
      TCSG=COS(PHI2*DEGRAD)
      DEN=X(4)*TSGN-X(5)*TCSG
      IF(DEN.EQ.0.0) GO TO 20
      SP=(X(2)*TCSG-X(1)*TSGN)/DEN
      IF(SP.LT.0.0) GO TO 20
      IF(ABS(TSGN).GT.1.E-6.AND.(X(2)+SP*X(5))*TSGN
     +       .LT.0.)GO TO 20
      IF(ABS(TCSG).GT.1.E-6.AND.(X(1)+SP*X(4))*TCSG
     +       .LT.0.)GO TO 20
C
      IF(SP.LT.SNXT) SNXT=SP
C
   20 CONTINUE
C
      TH2=PAR(4)
      IBOUN=0
   30 CONTINUE
      IBOUN=IBOUN+1
C
      STH=SIN(TH2*DEGRAD)
      CTH=COS(TH2*DEGRAD)
      STH2=STH*STH
      CTH2=CTH*CTH
      AA=STH2*X(6)**2-CTH2*(X(4)**2+X(5)**2)
      BB=STH2*X(6)*X(3)-CTH2*(X(4)*X(1)+
     +X(5)*X(2))
      CC=STH2*X(3)**2-CTH2*(X(1)**2+X(2)**2)
      SQ=BB*BB-AA*CC
      IF(SQ.LT.0.0) GO TO 60
      ST=BIG
      IF(ABS(BB).GE.SMALL) ST=-CC*0.5/BB
      ITRY=1
      IF(AA.EQ.0.0) GO TO 45
      T2=SQRT(SQ)/AA
      ITRY=0
      T1=-BB/AA
   40 CONTINUE
      ST=T1+T2
   45 CONTINUE
      IF((X(3)+ST*X(6))*CTH.LT.0.0) GO TO 50
      IF(ST.GT.0.0.AND.ST.LT.SNXT) SNXT=ST
   50 CONTINUE
      IF(ITRY.NE.0) GO TO 60
      ITRY=1
      T2=-T2
      GO TO 40
   60 CONTINUE
      IF(IBOUN.NE.1) GO TO 999
      TH2=PAR(3)
      GO TO 30
C
   70 CONTINUE
C              THIS BIT FOR X,Y,Z=0,0,0
C
      SNXT=PAR(2)
C              WE HAVE IGNORED THETA AND PHI BUT USERS SHOULDN'T
C              USE THETA PHI SEGMENTATION AT R=0 ANYWAY.
C
  999 CONTINUE
      END
