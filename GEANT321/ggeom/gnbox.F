*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNBOX (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'BOX' VOLUME,   *
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
C.    *         Authors  R.Brun,  P. Weidhaas  *********               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION X(6), PAR(6)
C.
C.    ------------------------------------------------------------------
C.
      SNXT = BIG
      S1   = BIG
      S2   = BIG
      S3   = BIG
      SAFX1= PAR(1)-X(1)
      SAFX2= PAR(1)+X(1)
      SAFY1= PAR(2)-X(2)
      SAFY2= PAR(2)+X(2)
      SAFZ1= PAR(3)-X(3)
      SAFZ2= PAR(3)+X(3)
 
      IF (IACT .LT. 3) THEN
 
C       -------------------------------------------------
C       |  Compute safety-distance 'SAFE' (P.Weidhaas)  |
C       -------------------------------------------------
 
        SAFE  = MIN(SAFX1,SAFX2, SAFY1,SAFY2, SAFZ1,SAFZ2)
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
      IF(X(4).GT.0.) S1 = SAFX1/X(4)
      IF(X(4).LT.0.) S1 =-SAFX2/X(4)
      IF(X(5).GT.0.) S2 = SAFY1/X(5)
      IF(X(5).LT.0.) S2 =-SAFY2/X(5)
      IF(X(6).GT.0.) S3 = SAFZ1/X(6)
      IF(X(6).LT.0.) S3 =-SAFZ2/X(6)
 
      SNXT  = MIN(S1,S2,S3)
 
  999 CONTINUE
      END
