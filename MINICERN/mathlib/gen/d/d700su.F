*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D700SU
C     PREPARE FIXED DATA FOR FAST FOURIER TRANSFORM.
C
      COMMON /D700DT/ N, N2, N4, M, F, RTTWO
      COMMON /FWORK/ W(321)
C
      IF(M .LT. 2) GO TO 90
      N4 = 2 ** (M-2)
      N2 = N4 + N4
      N  = N2 + N2
      F  = 1.0 / SQRT(REAL(N))
      DA = 4.0 * ATAN(1.0) / N2
      RTTWO = SQRT(2.0)
C
      NC = N4 - 1
      IF(NC .LE. 0) RETURN
      DO 10 MC = 1, NC
      W(MC) = COS(DA*MC)
   10 CONTINUE
C
      RETURN
C     M MUST BE AT LEAST 2, I.E. THE PERIOD HAS AT LEAST 4 POINTS.
   90 WRITE(6,100) M
      STOP
C
  100 FORMAT('1* ERROR IN RFT ... M =',I5,', SHOULD BE .GE. 2 *')
      END
