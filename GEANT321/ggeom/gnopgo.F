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
      SUBROUTINE GNOPGO (X, PAR, IACT, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       COMPUTE DISTANCE UP TO INTERSECTION WITH 'PGON' VOLUME,  *
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
C.    *                       from point X(1-3) along direction X(4-6) *
C.    *       SAFE  (output) : shortest distance from point X(1-3) to  *
C.    *                       any boundary (possibly overestimated)    *
C.    *                                                                *
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Authors   R.Nierhaus, F.Bruyant   *********            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      REAL X(*),PAR(*)
C
C.    ------------------------------------------------------------------
*
      SNXT = BIG
      IF (IACT.LT.3) THEN
*
* ***   Compute SAFE
*
         R2 = X(1)**2 +X(2)**2
         R  = SQRT(R2)
         DPHI2  = 0.5*PAR(2)*DEGRAD/PAR(3)
         CDPHI2 = COS(DPHI2)
         NZ = PAR(4)
*
         IF (X(3).LT.PAR(5)) THEN
            SAFEZ = PAR(5) -X(3)
         ELSE IF (X(3).GT.PAR(3*NZ+2)) THEN
            SAFEZ = X(3) -PAR(3*NZ+2)
         ELSE
            SAFEZ = BIG
         ENDIF
*
         RMIN = BIG
         RMAX = 0.
         DO 10 I = 6,3*NZ+3,3
            IF (PAR(I).LT.RMIN) RMIN = PAR(I)
            IF (PAR(I+1).GT.RMAX) RMAX = PAR(I+1)
   10    CONTINUE
         RMAX = RMAX/CDPHI2
*
*  **   Treat outside Z-limits first
*
         IF (SAFEZ.LT.0.9*BIG) THEN
            IF (R.LT.RMIN) THEN
               SAFE = MAX(SAFEZ,RMIN-R)
            ELSE IF (R.GT.RMAX) THEN
               SAFE = MAX(SAFEZ,R-RMAX)
            ELSE
               SAFE = SAFEZ
            ENDIF
            GO TO 50
         ELSE
*
*  **   Now one is inside the Z-limits
*
            DO 20 I = 8,3*NZ+2,3
              IF (X(3).GT.PAR(I)) GO TO 20
              IF (PAR(I+1).EQ.PAR(I-2)) THEN
                 IF (R.LT.PAR(I+1)) THEN
                    SAFE = PAR(I+1) -R
                    IF (NZ.EQ.2) GO TO 50
                    XT3  = X(3) -0.5*(PAR(I)+PAR(I-3))
                    DZ   = 0.5*(PAR(I)-PAR(I-3))
                    GO TO 18
                 ENDIF
              ELSE IF (PAR(I+2).EQ.PAR(I-1)) THEN
                 IF (R.GT.PAR(I+2)) THEN
                    SAFE = R -PAR(I+2)
                    IF (NZ.EQ.2) GO TO 50
                    XT3  = X(3) -0.5*(PAR(I)+PAR(I-3))
                    DZ   = 0.5*(PAR(I)-PAR(I-3))
                    GO TO 18
                 ENDIF
              ENDIF
*
              XT3  = X(3) -0.5*(PAR(I)+PAR(I-3))
              DZ   = 0.5*(PAR(I)-PAR(I-3))
              FACT = 0.5*(PAR(I+1)-PAR(I-2))/DZ
              RAD1 = PAR(I-2) +FACT*(DZ+XT3)
              IF (R.LT.RAD1) THEN
                 SAFE = (RAD1-R)/SQRT(1.+FACT**2)
              ELSE
                 FACT = 0.5*(PAR(I+2)-PAR(I-1))/DZ
                 RAD2 = PAR(I-1) +FACT*(DZ+XT3)
                 RAD2=RAD2/CDPHI2
                 IF (R.GT.RAD2) THEN
                    SAFE = (R-RAD2)/SQRT(1.+FACT**2)
                 ELSE
                    SAFE = 0.
                 ENDIF
              ENDIF
              IF (NZ.EQ.2) GO TO 19
*
   18         SAFEZ = DZ-ABS(XT3)
              SAFE  = MIN(SAFEZ,SAFE)
*
   19         IF (R.LT.RMIN) THEN
                 SAFE = MAX(RMIN-R,SAFE)
              ELSE IF (R.GT.RMAX) THEN
                 SAFE = MAX(R-RMAX,SAFE)
              ENDIF
              GO TO 50
   20       CONTINUE
         ENDIF
*
*  **   Case of a Phi segmented polygon
*               (not coded yet)
   50    CONTINUE
*
         IF (IACT .EQ. 0) GO TO 999
         IF (IACT .EQ. 1) THEN
            IF (SNEXT .LT. SAFE) GO TO 999
         ENDIF
      ENDIF
*
      CALL GNOPG1 (X, PAR, SNXT)
*
  999 END
