*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE GAMRN.FOR
*COPY GAMRN
*
*=== gamrn ============================================================*
*
      FUNCTION GAMRN(ALAM,ETA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      REAL RNDM(2)
C
C********************************************************************
C
C     RANDOM NUMBER SELECTION FROM GAMMA DISTRIBUTION
C     F(X) = ALAM**ETA*X**(ETA-1)*EXP(-ALAM*X) / GAM(ETA)
C
C********************************************************************
C
      NCOU = 0
      N = ETA
      F = ETA - N
      IF (F .EQ. 0.D0) GO TO 100
  10  CALL GRNDM(RNDM,1)
      R = RNDM(1)
      NCOU = NCOU + 1
      IF (NCOU.GE.11) GO TO 100
      IF (R .LT. F/(F + 2.718281828459045D0)) GO TO 20
      CALL GRNDM(RNDM,1)
      YYY = LOG(RNDM(1)+1.0D-38)/F
      IF (ABS(YYY) .GT. 50.D0) GO TO 100
      Y = EXP(YYY)
      CALL GRNDM(RNDM,1)
      IF (LOG(RNDM(1) + 1.0D-38) .GT. -Y) GO TO 10
      GO TO 50
 100  Y = 0.D0
      GO TO 70
  20  CALL GRNDM(RNDM,2)
      Y = 1.D0 - LOG(RNDM(1) + 1.0D-38)
      IF (RNDM(2) .GT. Y**(F-1.D0)) GO TO 10
  50  IF (N .EQ. 0) GO TO 150
  70  Z = 1.D0
      IF (N .LE. 0) N = 1
      DO 80 I = 1,N
  80  CALL GRNDM(RNDM,1)
      Z = Z * RNDM(1)
      Y = Y - LOG(Z + 1.0D-38)
 150  GAMRN = Y/ALAM
      RETURN
      END
