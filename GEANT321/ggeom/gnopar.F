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
      SUBROUTINE GNOPAR (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH 'PARA' VOLUME,   *
C.    *       FROM OUTSIDE POINT X(1-3) ALONG DISTANCE X(4-6)          *
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
      DIMENSION X(6), PAR(6)
C.
C.    -------------------------------------------
C.
      SNXT=BIG
      DN31=-PAR(3)-X(3)
      DN32=PAR(3)-X(3)
      YT=X(2)-PAR(6)*X(3)
      DN21=-PAR(2)-YT
      DN22=PAR(2)-YT
      CTY=1.0/SQRT(1.0+PAR(6)**2)
C
      XT=X(1)-PAR(5)*X(3)-PAR(4)*YT
      DN11=-PAR(1)-XT
      DN12=PAR(1)-XT
      CTX=1.0/SQRT(1.0+PAR(4)**2+PAR(5)**2)
C
      SN3=DN31
      IF(SN3.LT.0.0) SN3=-DN32
      SN2=DN21*CTY
      IF(SN2.LT.0.0) SN2=-DN22*CTY
      SN1=DN11*CTX
      IF(SN1.LT.0.0) SN1=-DN12*CTX
 
      IF (IACT .LT. 3) THEN
 
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (McPherson)   |
C       -------------------------------------------------
 
        SAFE=MAX(SN1,SN2,SN3)
C
C      IF(SAFETY.LT.0.0) GO TO 900
C
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
      IF(X(6).NE.0.0) GO TO 5
      IF(DN32*DN31.GT.0.0) GO TO 999
      DN31=0.0
      DN32=BIG
      GO TO 10
    5 CONTINUE
C
      DN31=DN31/X(6)
      DN32=DN32/X(6)
      IF(X(6).GT.0.0) GO TO 10
C
      DNS=DN31
      DN31=DN32
      DN32=DNS
C
   10 CONTINUE
C
      IF(DN32.LE.0.0) GO TO 999
C
      DY=X(5)-PAR(6)*X(6)
C
      IF(DY.NE.0.0) GO TO 15
      IF(DN21*DN22.GT.0.0) GO TO 999
      DN21=0.0
      DN22=BIG
      GO TO 20
   15 CONTINUE
      DN21=DN21/DY
      DN22=DN22/DY
      IF(DY.GT.0.0) GO TO 20
C
      DNS=DN21
      DN21=DN22
      DN22=DNS
C
   20 CONTINUE
C
      IF(DN22.LE.0.0) GO TO 999
C
      DX=X(4)-PAR(5)*X(6)-PAR(4)*DY
      IF(DX.NE.0.0) GO TO 25
      IF(DN11*DN12.GT.0.0) GO TO 999
      DN11=0.0
      DN12=BIG
      GO TO 30
   25 CONTINUE
      DN11=DN11/DX
      DN12=DN12/DX
      IF(DX.GT.0.0) GO TO 30
C
      DNS=DN11
      DN11=DN12
      DN12=DNS
C
   30 CONTINUE
C
      IF(DN12.LE.0.0) GO TO 999
C
      SMIN=MAX(DN11,DN21,DN31)
      SMAX=MIN(DN12,DN22,DN32)
C
      IF(SMAX.LE.SMIN) GO TO 999
C
      IF(SMIN.LE.0.)GO TO 999
      SNXT = SMIN
C
  999 CONTINUE
      END
