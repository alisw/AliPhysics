*
* $Id$
*
* $Log$
* Revision 1.2  1997/07/01 15:31:41  mclareni
* Correction for FP exceptions and restriction of calculations to realistic values, from M. Schroder
*
* Revision 1.1.1.1  1996/04/01 15:01:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION ERF(X)

      LOGICAL LEF
      DIMENSION P2(0:4),Q2(0:4)

      PARAMETER(Z1 = 1, HF = Z1/2, C1 = 0.56418 958)

C     Above the value of VMAX any calculation is pointless. The value is
C     choosen with a big safety margin - even the double precision
C     version only returns 1. for V (=ABS(X)) >= 5.9
      PARAMETER(VMAX = 7.)
C     The value for SWITCH is badly chosen for the single precision
C     version, which returns 1. already for V >= 3.9
      PARAMETER(SWITCH = 4.)

      DATA P10,Q10,P11 /+3.67678 77, +3.25845 93, -9.79704 65E-2/

      DATA (P2(I),Q2(I),I=0,4)
     +/+7.37388 83E+0, +7.37396 09E+0, +6.86501 85E+0, +1.51849 08E+1,
     1 +3.03179 93E+0, +1.27955 30E+1, +5.63169 62E-1, +5.35421 68E+0,
     2 +4.31877 87E-5, +1/

      DATA P30,Q30,P31 /-1.24368 54E-1, +4.40917 06E-1, -9.68210 36E-2/

      LEF=.TRUE.
      GO TO 9

      ENTRY ERFC(X)
      LEF=.FALSE.

    9 V=ABS(X)
      IF(V .LT. HF) THEN
       Y=V**2
       H=X*(P10+P11*Y)/(Q10+Y)
       HC=1-H
      ELSE
       IF(V .LT. SWITCH) THEN
        AP=P2(4)
        AQ=Q2(4)
        DO 2 I = 3,0,-1
        AP=P2(I)+V*AP
    2   AQ=Q2(I)+V*AQ
        HC=EXP(-V**2)*AP/AQ
        H=1-HC
       ELSEIF ( V .LT. VMAX) THEN
        Y=1/V**2
        HC=EXP(-V**2)*(C1+Y*(P30+P31*Y)/(Q30+Y))/V
        H=1-HC
C     for very big values we can save us any calculation, and the
C     FP-exceptions we would get from EXP.
       ELSE
        H  = 1.
        HC = 0.
       ENDIF
       IF(X .LT. 0) THEN
        H=-H
        HC=2-HC
       ENDIF
      ENDIF
      IF(LEF) THEN
       ERF=H
      ELSE
       ERFC=HC
      ENDIF
      RETURN
      END
#endif
