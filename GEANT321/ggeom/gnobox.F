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
      SUBROUTINE GNOBOX (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'BOX' VOLUME,   *
C.    *        FROM OUTSIDE POINT X(1-3) ALONG DIRECTION X(4-6)        *
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
C.    *       ==>Called by : GNEXT, GTNEXT                             *
C.    *       Author  A.McPherson, P.Weidhaas  ******                  *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6),PAR(3),SAF(3)
      EQUIVALENCE (SAF(1),SAF1),(SAF(2),SAF2),(SAF(3),SAF3)
C.
C.    ----------------------------------------------------------------
C.
 
      SNXT = BIG
      DO 5 I=1,3
        SAF(I)=ABS(X(I))-PAR(I)
   5  CONTINUE
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
        SAFE  = MAX (SAF1, SAF2, SAF3, 0.0)
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
      IF(X(1)*X(4).GE.0.0.AND.SAF1.GT.0.) GO TO 999
      IF(X(2)*X(5).GE.0.0.AND.SAF2.GT.0.) GO TO 999
      IF(X(3)*X(6).GE.0.0.AND.SAF3.GT.0.) GO TO 999
C
      SMAX1=BIG
      SMAX2=BIG
      SMAX3=BIG
      SMIN1=0.0
C
      IF(SAF1.LT.0.) GO TO 10
      IF(X(4).EQ.0.0) GO TO 999
      S4 = ABS(X(4))
      SMIN1=SAF1 / S4
      SMAX1=(PAR(1)+ABS(X(1))) / S4
      GO TO 20
C
   10 CONTINUE
      IF(X(4).NE.0.)THEN
         SMAX1=PAR(1)/ABS(X(4))-X(1)/X(4)
      ENDIF
   20 CONTINUE
C
      SMIN2=0.0
      IF(SAF2.LT.0.) GO TO 30
      IF(X(5).EQ.0.0) GO TO 999
      S5=ABS(X(5))
      SMIN2=SAF2/S5
      SMAX2=(PAR(2)+ABS(X(2)))/S5
      GO TO 40
C
   30 CONTINUE
      IF(X(5).NE.0.)THEN
         SMAX2=PAR(2)/ABS(X(5))-X(2)/X(5)
      ENDIF
   40 CONTINUE
C
      SMIN3=0.0
      IF(SAF3.LT.0.) GO TO 50
      IF(X(6).EQ.0.0) GO TO 999
      S6=ABS(X(6))
      SMIN3=SAF3/S6
      SMAX3=(PAR(3)+ABS(X(3)))/S6
      GO TO 60
C
   50 CONTINUE
      IF(X(6).NE.0.)THEN
         SMAX3=PAR(3)/ABS(X(6))-X(3)/X(6)
      ENDIF
   60 CONTINUE
C
      SMIN=MAX(SMIN1,SMIN2,SMIN3)
      SMAX=MIN(SMAX1,SMAX2,SMAX3)
      IF(SMAX.LT.SMIN) GO TO 999
C
      IF(SMIN.LE.0.)GO TO 999
      SNXT = SMIN
 
  999 CONTINUE
      END
