*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION  ATANI(X)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DATANI(X)
#include "gen/imp64.inc"
#endif
      DIMENSION C(0:18)

      PARAMETER (PI = 3.14159 26535 89793 24D0, PIH = PI/2)

      DATA C( 0) / 0.93432 00492 92895 95D0/
      DATA C( 1) /-0.01950 47944 34351 90D0/
      DATA C( 2) / 0.00125 17037 06300 59D0/
      DATA C( 3) /-0.00011 19241 41205 64D0/
      DATA C( 4) / 0.00001 17754 53855 07D0/
      DATA C( 5) /-0.00000 13652 83304 87D0/
      DATA C( 6) / 0.00000 01688 82892 74D0/
      DATA C( 7) /-0.00000 00218 80246 51D0/
      DATA C( 8) / 0.00000 00029 35063 61D0/
      DATA C( 9) /-0.00000 00004 04523 41D0/
      DATA C(10) / 0.00000 00000 56976 00D0/
      DATA C(11) /-0.00000 00000 08168 83D0/
      DATA C(12) / 0.00000 00000 01188 69D0/
      DATA C(13) /-0.00000 00000 00175 16D0/
      DATA C(14) / 0.00000 00000 00026 09D0/
      DATA C(15) /-0.00000 00000 00003 92D0/
      DATA C(16) / 0.00000 00000 00000 59D0/
      DATA C(17) /-0.00000 00000 00000 09D0/
      DATA C(18) / 0.00000 00000 00000 01D0/

      V=ABS(X)
      IF(V .LE. 1) THEN
       R=0
       Y=V
      ELSE
       R=PIH*LOG(V)
       Y=1/V
      END IF
      H=2*Y**2-1
      ALFA=H+H
      B1=0
      B2=0
      DO 1 I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      R=R+Y*(B0-B2)
      IF(X .LT. 0) R=-R
#if !defined(CERNLIB_DOUBLE)
       ATANI=R
#endif
#if defined(CERNLIB_DOUBLE)
      DATANI=R
#endif
      RETURN
      END
