*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:55  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CGAMMA(Z)
#include "gen/defc64.inc"
     + CGAMMA
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION WGAMMA(Z)
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +  WGAMMA
#endif
#include "gen/defc64.inc"
     +       Z,U,V,F,H,S
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CGAMMA')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CGAMMA/WGAMMA')
#endif
      DIMENSION C(0:15)

      PARAMETER (Z1 = 1, HF = Z1/2)

#if defined(CERNLIB_QF2C)
#include "gen/gcmpfun.inc"
#endif

      DATA PI /3.14159 26535 89793 24D0/
      DATA C1 /2.50662 82746 31000 50D0/

      DATA C( 0) / 41.62443 69164 39068D0/
      DATA C( 1) /-51.22424 10223 74774D0/
      DATA C( 2) / 11.33875 58134 88977D0/
      DATA C( 3) / -0.74773 26877 72388D0/
      DATA C( 4) /  0.00878 28774 93061D0/
      DATA C( 5) / -0.00000 18990 30264D0/
      DATA C( 6) /  0.00000 00019 46335D0/
      DATA C( 7) / -0.00000 00001 99345D0/
      DATA C( 8) /  0.00000 00000 08433D0/
      DATA C( 9) /  0.00000 00000 01486D0/
      DATA C(10) / -0.00000 00000 00806D0/
      DATA C(11) /  0.00000 00000 00293D0/
      DATA C(12) / -0.00000 00000 00102D0/
      DATA C(13) /  0.00000 00000 00037D0/
      DATA C(14) / -0.00000 00000 00014D0/
      DATA C(15) /  0.00000 00000 00006D0/

#if !defined(CERNLIB_QF2C)
#include "gen/gcmpfun.inc"
#endif

      U=Z
      X=U
      IF(GIMAG(U) .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       F=0
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C305.1',ERRTXT)
      ELSE
       IF(X .GE. 1) THEN
        F=1
        V=U
       ELSEIF(X .GE. 0) THEN
        F=1/U
        V=1+U
       ELSE
        F=1
        V=1-U
       END IF
       H=1
       S=C(0)
       DO 1 K = 1,15
       H=((V-K)/(V+(K-1)))*H
    1  S=S+C(K)*H
       H=V+(4+HF)
       H=C1*EXP((V-HF)*LOG(H)-H)*S
       IF(X .LT. 0) H=PI/(SIN(PI*U)*H)
      ENDIF
#if !defined(CERNLIB_DOUBLE)
      CGAMMA=F*H
#endif
#if defined(CERNLIB_DOUBLE)
      WGAMMA=F*H
#endif
      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
