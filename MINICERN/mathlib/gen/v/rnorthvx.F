*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_VAX)
      REAL*4 FUNCTION RNORTH(K)
C
C RNOR TOOTH FUNCTION
C
      IMPLICIT        NONE
C
      INTEGER*4       K
C
      REAL*4          C (45)
C
      DATA C/'2B5F407D'X, '2B5F407D'X, 'A9AD407A'X, 'A6484075'X,
     1       '24964073'X, '2131406E'X, '9C1A4066'X, '98B54061'X,
     2       '139E405A'X, '8E874052'X, '87BE4048'X, '02A64041'X,
     3       'FBDD4036'X, 'F513402C'X, 'EE4A4022'X, 'E7804018'X,
     4       '62694011'X, '5BA04007'X, 'A9AC3FFA'X, '9C1A3FE6'X,
     5       '91EC3FD7'X, '84583FC3'X, '7A2A3FB4'X, '6FFC3FA5'X,
     6       '65CE3F96'X, '5BA03F87'X, 'A2E43F70'X, '95503F5C'X,
     7       '80F43F3E'X, '73603F2A'X, '65CC3F16'X, '58383F02'X,
     8       'A2E03EF0'X, '87B83EC8'X, '7A283EB4'X, '6C983EA0'X,
     9       '5F003E8C'X, 'A2E43E70'X, '87BE3E48'X, '6C983E20'X,
     1       'A2E43DF0'X, 'A2E43DF0'X, '6C983DA0'X, '6C983DA0'X,
     2       '6C983DA0'X/
C
      INTEGER*4       I1      /'FBC35400'X/
      INTEGER*4       I2      /'FE79702E'X/
C
      REAL*4          S, T, B, X, P, V
      INTEGER*4       J
C
      REAL*4 VNI, UNI
      EXTERNAL        VNI, UNI
C
      IF (K .GT. I1) GOTO 3
      S = UNI (0)
      T = UNI (0)
      B = AINT (7. * (S + T) + 37. * ABS (S - T))
      X = UNI (0) - UNI (0)
      RNORTH = .0625 * (X + SIGN (B,X))
      RETURN
 3    IF (K .GT. I2) GOTO 5
 4    RNORTH = 2.75 * VNI (0)
      J = 16. * ABS (RNORTH) + 1.
      IF (J - 14) 6,6,7
 6    P = (J + J - 1) * .1497466E-2
      GOTO 8
 7    P= (89 - J - J) * .698817E-3
 8    IF (UNI (0) .GT. 79.78846 * (EXP (-.5 * RNORTH * RNORTH)
     1       -C (J) - P * (J - 16. * ABS (RNORTH)))) GOTO 4
      RETURN
 5    V = VNI (0)
      IF (V .EQ. 0) GOTO 5
      X = SQRT (7.5625 - 2. * LOG (ABS (V)))
      IF (UNI (0) * X .GT. 2.75) GOTO 5
      RNORTH = SIGN (X, V)
      RETURN
      END
#endif
