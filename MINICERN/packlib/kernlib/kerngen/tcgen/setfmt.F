*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_TCGEN)
      SUBROUTINE SETFMT (FMTLET,FMTNUM,NDIG,XX,NX)
C
C CERN PROGLIB# M224    SETFMT          .VERSION KERNFOR  4.08  840613
C ORIG. 01/08/71
C
      DIMENSION    XX(9)
      INTEGER      FMTLET,FMTNUM, FMTN(10),EE,EF
      DATA  FMTN   / 4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,4H5   ,
     +                      4H6   ,4H7   ,4H8   ,4H9    /
      DATA  EE,EF  / 4HE   ,4HF    /
      DATA  VERYSM / 1.E-36 /
C
C
      XBIG= 0.
      MXV = MIN (NDIG,9)
      N   = 0
C
      DO 9 J=1,NX
    9 XBIG =  MAX (ABS(XX(J)), XBIG)
C
      IF (XBIG .EQ. 0.)      GO TO 24
      IF (XBIG .LT. VERYSM)  GO TO 21
C
C----      NINT = NO. OF DIGITS BEFORE THE DECIMAL POINT
C--        N    = NO. OF DIGITS AFTER  THE DECIMAL POINT
C--       -NINT = NO.OF ZEROES AFTER THE DECIMAL POINT, IF PURE FRACTION
C
      NINT = INT (LOG10(XBIG)+100.) - 99
      IF  (NINT .GT. NDIG)  GO TO 21
      IF (-NINT .GT. MXV-2) GO TO 21
      MXV= MIN(9, MAX(0, NDIG-MAX(0,NINT)) )
      BIAS= .25*10.**(-MXV)
      TOL = BIAS+BIAS
C
      DO 19 J=1,NX
      X = (ABS(XX(J))+BIAS) * 10.**N
C
   12 IF (N .GE. MXV)  GO TO 24
      X = X - AINT(X)
      IF (X .LT. TOL)  GO TO 19
      TOL= 10.*TOL
      X  = 10.*X
      N  = N+1
      GO TO 12
   19 CONTINUE
      GO TO 24
C
   21 FMTLET= EE
      N     = MIN(5, NDIG-5)
      GO TO 25
C
   24 FMTLET= EF
   25 FMTNUM= FMTN(N+1)
      RETURN
C
      END
#endif
