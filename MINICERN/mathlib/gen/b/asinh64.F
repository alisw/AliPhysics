*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DASINH(X)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='DASINH')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION ASINH(X)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='ASINH')
#endif
      DIMENSION C(0:19)

      DATA C( 0) / 0.90649 39198 46333 18D0/
      DATA C( 1) /-0.02704 21478 78869 64D0/
      DATA C( 2) / 0.00211 68145 57973 56D0/
      DATA C( 3) /-0.00021 76650 54603 40D0/
      DATA C( 4) / 0.00002 55196 04364 81D0/
      DATA C( 5) /-0.00000 32329 14485 29D0/
      DATA C( 6) / 0.00000 04310 66959 88D0/
      DATA C( 7) /-0.00000 00596 06134 55D0/
      DATA C( 8) / 0.00000 00084 69211 32D0/
      DATA C( 9) /-0.00000 00012 29008 59D0/
      DATA C(10) / 0.00000 00001 81376 79D0/
      DATA C(11) /-0.00000 00000 27138 46D0/
      DATA C(12) / 0.00000 00000 04107 37D0/
      DATA C(13) /-0.00000 00000 00627 70D0/
      DATA C(14) / 0.00000 00000 00096 72D0/
      DATA C(15) /-0.00000 00000 00015 01D0/
      DATA C(16) / 0.00000 00000 00002 34D0/
      DATA C(17) /-0.00000 00000 00000 37D0/
      DATA C(18) / 0.00000 00000 00000 06D0/
      DATA C(19) /-0.00000 00000 00000 01D0/

      V=ABS(X)
      IF(V .LE. 1) THEN
       H=2*V**2-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       R=SIGN(V*(B0-B2),X)
      ELSE
       R=LOG(X+SQRT(1+X**2))
      ENDIF
#if defined(CERNLIB_DOUBLE)
      DASINH=R
#endif
#if !defined(CERNLIB_DOUBLE)
       ASINH=R
#endif
      RETURN
      END
