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
      SUBROUTINE GNTRAP (X, PAR, IACT, IFLAG, SNEXT, SNXT, SAFE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMPUTE DISTANCE UP TO INTERSECTION WITH 'TRD1' OR 'TRD2' *
C.    *      VOLUME, FROM INSIDE POINT X(1-3) ALONG DIRECTION X(4-6)   *
C.    *                                                                *
C.    *       PAR   (input)  : volume parameters                       *
C.    *       IACT  (input)  : action flag                             *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, and SNXT only if SNEXT .GT.new SAFE *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *       IFLAG (input)  : 1 for TRD1, 2 for TRD2                  *
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
C.    ------------------------------------------------------------------
C.
      SNXT = BIG
      IF (IFLAG .EQ. 1) THEN
        FACTX  = (PAR(2) - PAR(1)) / (2.0 * PAR(4))
      ELSEIF (IFLAG .EQ. 2) THEN
        FACTX  = (PAR(2) - PAR(1)) / (2.0 * PAR(5))
        FACTY  = (PAR(4) - PAR(3)) / (2.0 * PAR(5))
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
 
          SAF2  = PAR(3) - ABS(X(2))
          SAF3  = PAR(4) - ABS(X(3))
C
C......  Distance alng x-direction to slanted wall.
C
          DISTX = PAR(1) + FACTX * (PAR(4) + X(3)) - ABS(X(1))
          SAF1  = DISTX / SQRT (1.0 + FACTX*FACTX)
        ELSE
C*******************************************************
C
C......  Shape "TRD2":   both x and y vary with z.
C
C*******************************************************
          SAF3  = PAR(5) - ABS(X(3))
C
C......  Distance along x-direction to slanted wall.
C
          DISTX = PAR(1) + FACTX * (PAR(5) + X(3)) - ABS(X(1))
          SAF1  = DISTX / SQRT (1.0 + FACTX*FACTX)
C
C......  Distance along y-direction to slanted wall.
C
          DISTY = PAR(3) + FACTY * (PAR(5) + X(3)) - ABS(X(2))
          SAF2  = DISTY / SQRT (1.0 + FACTY*FACTY)
        ENDIF
 
        SAFE  = MIN (SAF1, SAF2, SAF3)
 
        IF (IACT .EQ. 0) GO TO 999
        IF (IACT .EQ. 1) THEN
          IF (SNEXT .LT. SAFE) GO TO 999
        ENDIF
      ENDIF
 
 
C     ------------------------------------------------
C     |  Compute vector-distance 'SNXT' (McPherson)  |
C     ------------------------------------------------
 
 
      IF(IFLAG.NE.1) GO TO 20
C
C              Only x-thickness varies with z.
C
C
C              First check z.
C
      ZEND=PAR(4)
      IF(X(6).LT.0.0) ZEND=-ZEND
      IF(X(6).EQ.0.)GO TO 5
      SNXT=(ZEND-X(3))/X(6)
C
C              Now do Y.
C
   5  YEND=PAR(3)
      IF(X(5).LT.0.0) YEND=-YEND
      IF(X(5).EQ.0.)GO TO 7
      SN=(YEND-X(2))/X(5)
      IF(SN.LT.SNXT) SNXT=SN
C
C              Now do X.
C
   7  DXM=0.5*(PAR(1)+PAR(2))
      ANUM=DXM+FACTX*X(3)-X(1)
      DENO=X(4)-FACTX*X(6)
C
      IF(DENO.EQ.0.0) GO TO 10
      QUOT = ANUM/DENO
      IF(QUOT.LE.0.0.OR.QUOT.GT.SNXT) GO TO 10
      SNXT=QUOT
      GO TO 999
C
   10 CONTINUE
C
      ANUM=-FACTX*X(3)-X(1)-DXM
      DENO=FACTX*X(6)+X(4)
C
      IF(DENO.EQ.0.0) GO TO 999
      QUOT = ANUM/DENO
      IF(QUOT.GT.0.0.AND.QUOT.LT.SNXT) SNXT=QUOT
C
      GO TO 999
   20 CONTINUE
C
C              x- and y-thicknesses vary with z.
C
C
C              First check z
C
      ZEND=PAR(5)
      IF(X(6).LT.0.0) ZEND=-ZEND
      IF(X(6).EQ.0.)GO TO 25
      SNXT=(ZEND-X(3))/X(6)
C
C              Now do x.
C
  25  DXM=0.5*(PAR(1)+PAR(2))
      ANUM=DXM+FACTX*X(3)-X(1)
      DENO=X(4)-FACTX*X(6)
C
      IF(DENO.EQ.0.0) GO TO 40
      QUOT = ANUM/DENO
      IF(QUOT.LE.0.0.OR.QUOT.GT.SNXT) GO TO 40
      SNXT=QUOT
      GO TO 60
C
   40 CONTINUE
C
      ANUM=-FACTX*X(3)-X(1)-DXM
      DENO=FACTX*X(6)+X(4)
C
      IF(DENO.EQ.0.0) GO TO 60
      QUOT = ANUM/DENO
      IF(QUOT.GT.0.0.AND.QUOT.LT.SNXT) SNXT=QUOT
C
   60 CONTINUE
C
C             Now do y.
C
      DYM=0.5*(PAR(3)+PAR(4))
      ANUM=DYM+FACTY*X(3)-X(2)
      DENO=X(5)-FACTY*X(6)
C
      IF(DENO.EQ.0.0) GO TO 80
      QUOT = ANUM/DENO
      IF(QUOT.LE.0.0.OR.QUOT.GT.SNXT) GO TO 80
      SNXT=QUOT
      GO TO 999
C
   80 CONTINUE
C
      ANUM=-FACTY*X(3)-X(2)-DYM
      DENO=FACTY*X(6)+X(5)
C
      IF(DENO.EQ.0.0) GO TO 999
      QUOT = ANUM/DENO
      IF(QUOT.GT.0.0.AND.QUOT.LT.SNXT) SNXT=QUOT
C
  999 CONTINUE
      END
