*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNTRP (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'TRAP' VOLUME,  *
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
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6), PAR(11)
#if defined(CERNLIB_SINGLE)
      DIMENSION XT(6), P(11)
#endif
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION XT(6), P(11)
      DOUBLE PRECISION A,B,C,DSN,TOP1,TOP2,BOT1,BOT2
      DOUBLE PRECISION DX1,DX2,X1,X2,X3,DCX1,DCX2,DX,DY
      DOUBLE PRECISION SP1,SP2,DCX,DCY,TANG1,TANG2,A1,B1,A2,B2
      DOUBLE PRECISION XL1,XL2,XL3,XL4,CY,SY,SX,CX,C1,C2
      DOUBLE PRECISION TX,H0,DHDZ,CH
      DOUBLE PRECISION CHTAN,SHTAN
      DOUBLE PRECISION SN1,SN2,SN3,SN4,SN5,SN6,SN7,SN8,P15
#endif
C.
C.               --------------------------------------
C.
      SNXT = BIG
      SAFE = 0.0
      DO 1 I=1,6
   1  XT(I)=X(I)
      DO 2 I=1,11
   2  P(I)=PAR(I)
 
      P15=0.5/P(1)
      DHDZ=(P(4)-P(8))*P15
 
      IF (IACT .LT. 3) THEN
 
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE'  (McPherson)  |
C       -------------------------------------------------
 
 
C            CALCULATE RECTANGLE ON FACE AT Z=-DZ.
 
        DX1=P(4)*P(7)
        X1=-DX1-P(5)
        X2=DX1-P(6)
        IF(X2.GT.X1) X1=X2
        X2=-DX1+P(5)
        X3=DX1+P(6)
        IF(X3.LT.X2) X2=X3
        DCX1=(X1+X2)*0.5
        DX1=(X2-X1)*0.5
C
C            CALCULATE RECTANGLE ON FACE AT Z=DZ.
C
        DX2=P(8)*P(11)
        X1=-DX2-P(9)
        X2=DX2-P(10)
        IF(X2.GT.X1) X1=X2
        X2=-DX2+P(9)
        X3=DX2+P(10)
        IF(X3.LT.X2) X2=X3
        DCX2=(X1+X2)*0.5
        DX2=(X2-X1)*0.5
 
C
        TX=P(2)+(DCX2-DCX1)*P15
C
C            CALCULATE LOCAL RECTANGLE.
C
        SP1=(P(1)-XT(3))*P15
        SP2=(P(1)+XT(3))*P15
C
        DY=P(4)*SP1+P(8)*SP2
        DX=DX1*SP1+DX2*SP2
        DCX=(DCX1+DCX2)*0.5+TX*XT(3)
        DCY=P(3)*XT(3)
C
C            CHECK POINT IS INSIDE RECTANGLE.
C
        IF(ABS(XT(1)-DCX).GT.DX) GO TO 10
C
C            CALCULATE ANGLE OF YZ PLANES.
C
        TANG1=TX+(DX1-DX2)*P15
        TANG2=TX-(DX1-DX2)*P15
C
C            CALCULATE SAFETY FROM YZ PLANES.
C
        SAF1=(XT(1)-DCX+DX)/SQRT(1.0+TANG1*TANG1)
        SAF2=(DX-XT(1)+DCX)/SQRT(1.0+TANG2*TANG2)
C
C            CALCULATE ANGLE OF XZ PLANES.
C
        TANG1=P(3)+DHDZ
        TANG2=P(3)-DHDZ
C
C            CALCULATE SAFETY FROM XZ PLANES.
C
        SAF3=(XT(2)-DCY+DY)/SQRT(1.0+TANG1*TANG1)
        SAF4=(DY-XT(2)+DCY)/SQRT(1.0+TANG2*TANG2)
C
C            CALCULATE SAFETY FROM XY PLANES.
C
        SAF5=P(1)-ABS(XT(3))
C
C            OVERALL SAFETY.
C
        SAFE = MIN(SAF1,SAF2,SAF3,SAF4,SAF5)
 
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
   10 CONTINUE
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
 
C           FIRST FIND S TO Z LIMITS.
 
      SN1=BIG
      SN2=BIG
      SN3=BIG
      SN4=BIG
C
      IF(XT(6).EQ.0.0) GOTO 15
      SN1=(P(1)-XT(3))/XT(6)
      SN2=-(P(1)+XT(3))/XT(6)
      IF(SN1.GE.0.) GOTO 15
      SNX=SN1
      SN1=SN2
      SN2=SNX
15    CONTINUE
C
C           NOW Y LIMITS.
C
      H0=(P(4)+P(8))*0.5
C
      TOP1=H0-DHDZ*XT(3)
      TOP2=XT(2)-XT(3)*P(3)
      BOT1=XT(5)-XT(6)*P(3)
      BOT2=DHDZ*XT(6)
C
      IF(BOT1+BOT2.NE.0.0) SN3=(TOP1-TOP2)/(BOT1+BOT2)
      IF(BOT2-BOT1.NE.0.0) SN4=(TOP1+TOP2)/(BOT2-BOT1)
C
C           NOW X LIMTS.
C
      XL1=(P(5)+P(6)+P(9)+P(10))*0.25
      XL2=(P(5)-P(6)+P(9)-P(10))*0.25
      XL3=(P(5)+P(6)-P(9)-P(10))*0.25
      XL4=(P(5)-P(6)-P(9)+P(10))*0.25
C
      CY=XT(2)-XT(3)*P(3)
      SY=XT(5)-XT(6)*P(3)
      CH=H0-DHDZ*XT(3)
C
      A1=XL4*XT(6)*SY+XL3*XT(6)*XT(6)*DHDZ
      B1=XL4*(XT(6)*CY+XT(3)*SY)-XL3*XT(6)*(CH-XT(3)*DHDZ)-XL2*P(1)*SY-
     +XL1*XT(6)*DHDZ*P(1)
      C1=XL4*CY*XT(3)-XL3*XT(3)*CH-XL2*CY*P(1)+XL1*CH*P(1)
C
      CHTAN=(P(7)*P(4)+P(11)*P(8))*0.5
      SHTAN=-(P(7)*P(4)-P(11)*P(8))*P15
      CHTAN=CHTAN+SHTAN*XT(3)
      SHTAN=SHTAN*XT(6)
      CX=XT(1)-XT(3)*P(2)
      SX=XT(4)-XT(6)*P(2)
C
      A2=-P(1)*(DHDZ*XT(6)*SX+SY*SHTAN)
      B2=P(1)*(CH*SX-DHDZ*XT(6)*CX-CY*SHTAN-SY*CHTAN)
      C2=P(1)*(CH*CX-CY*CHTAN)
C
      ISN56=0
      A=A1+A2
      B=B1+B2
      C=C1+C2
      IF(B*B-A*C*4.0.LT.0.0) GO TO 20
      IF(ABS(A).LT.1.0E-7) GO TO 19
      DSN=SQRT(B*B-A*C*4.0)
      SN5=(-B+DSN)*0.5/A
      SN6=(-B-DSN)*0.5/A
      ISN56=1
      GO TO 20
   19 CONTINUE
      IF(ABS(B).LT.1.0E-5) GO TO 20
      ISN56=1
      SN5=-C/B
      SN6=-C/B
   20 CONTINUE
C
      ISN78=0
      A=A1-A2
      B=B1-B2
      C=C1-C2
      IF(B*B-A*C*4.0.LT.0.0) GO TO 30
      DSN=SQRT(B*B-A*C*4.0)
      IF(ABS(A).LT.1.0E-7) GO TO 25
      SN7=(-B+DSN)*0.5/A
      SN8=(-B-DSN)*0.5/A
      ISN78=1
      GO TO 30
   25 CONTINUE
      IF(ABS(B).LT.1.0E-5) GO TO 30
      ISN78=1
      SN7=-C/B
      SN8=-C/B
   30 CONTINUE
C
      IF(SN2.GT.0.0.AND.SN2.LT.SN1) SN1=SN2
      IF(SN3.GT.0.0.AND.SN3.LT.SN1) SN1=SN3
      IF(SN4.GT.0.0.AND.SN4.LT.SN1) SN1=SN4
      IF(ISN56.EQ.0) GO TO 40
      IF(SN5.GT.0.0.AND.SN5.LT.SN1) SN1=SN5
      IF(SN6.GT.0.0.AND.SN6.LT.SN1) SN1=SN6
   40 CONTINUE
      IF(ISN78.EQ.0) GO TO 50
      IF(SN7.GT.0.0.AND.SN7.LT.SN1) SN1=SN7
      IF(SN8.GT.0.0.AND.SN8.LT.SN1) SN1=SN8
   50 CONTINUE
C
      SNXT=SN1
C
  999 CONTINUE
      END
