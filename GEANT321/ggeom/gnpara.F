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
      SUBROUTINE GNPARA (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'PARA' VOLUME,  *
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
      DIMENSION X(6), PAR(6)
C.
C.    ---------------------------------------------------
C.
      SNXT = BIG
      DN31=PAR(3)-X(3)
      DN32=-PAR(3)-X(3)
C
      YT=X(2)-PAR(6)*X(3)
      DN21=PAR(2)-YT
      DN22=-PAR(2)-YT
      CTY=1.0/SQRT(1.0+PAR(6)**2)
C
      XT=X(1)-PAR(5)*X(3)-PAR(4)*YT
      DN11=PAR(1)-XT
      DN12=-PAR(1)-XT
      CTX=1.0/SQRT(1.0+PAR(4)**2+PAR(5)**2)
 
      IF (IACT .LT. 3) THEN
 
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (McPherson)   |
C       -------------------------------------------------
 
 
        SAFE=MIN(DN11*CTX,-DN12*CTX,DN21*CTY,-DN22*CTY,DN31,-DN32)
C
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
      SN3=BIG
      SN2=BIG
      SN1=BIG
C
      IF(X(6).NE.0.0) SN3=DN31/X(6)
      IF(SN3.LT.0.0) SN3=DN32/X(6)
C
      DY=X(5)-PAR(6)*X(6)
      IF(DY.NE.0.0) SN2=DN21/DY
      IF(SN2.LT.0.0) SN2=DN22/DY
C
      DX=X(4)-PAR(5)*X(6)-PAR(4)*DY
      IF(DX.NE.0.0) SN1=DN11/DX
      IF(SN1.LT.0.0) SN1=DN12/DX
C
      SNXT=MIN(SN1,SN2,SN3)
C
  999 CONTINUE
      END
