*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
C     This corresponds to CHISIN,IF=DOUBLE and CHISIN64,IF=-DOUBLE
      FUNCTION CHISIN(Q,N)

C     Computes the inverse of the
C     chi**2 distribution with n degrees of freedom.
C     Based on R.B. Goldstein,  ALGORITHM 451 Chi-Square Quantiles,
C     Collected Algorithms from CACM
C     Note the complementary definition of the integral!

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'CHISIN')
      DIMENSION C1(7),C2(7),C3(7),A(19)
      PARAMETER (Z1 = 1, HF = Z1/2)

      DATA C1

     1/ 1.565326E-3, 1.060438E-3,-6.950356E-3,-1.323293E-2, 2.277679E-2,
     2 -8.986007E-3,-1.513904E-2/
      DATA C2
     1/-1.450117E-3, 2.530010E-3, 5.169654E-3,-1.153761E-2, 1.128186E-2,
     2  2.607083E-2,-2.237368E-1/
      DATA C3
     1/ 9.780499E-5,-8.426812E-4, 3.125580E-3,-8.553069E-3, 1.348028E-4,
     2  4.713941E-1, 1.0000886/

      DATA A
     1/ 1.264616E-2,-1.425296E-2, 1.400483E-2,-5.886090E-3,-1.091214E-2,
     2 -2.304527E-2, 3.135411E-3,-2.728484E-4,-9.699681E-3, 1.316872E-2,
     3  2.618914E-2,-2.222222E-1, 5.406674E-5, 3.483789E-5,-7.274761E-4,
     4  3.292181E-3,-8.729713E-3, 4.714045E-1, 1/

      IF(Q .LT. 0 .OR. Q .GE. 1) THEN
       H=0
       WRITE(ERRTXT,101) Q
       CALL MTLPRT(NAME,'G101.1',ERRTXT)
      ELSEIF(N .LE. 0) THEN
       H=0
       WRITE(ERRTXT,102) N
       CALL MTLPRT(NAME,'G101.2',ERRTXT)
      ELSEIF(Q .EQ. 0) THEN
       H=0
      ELSEIF(N .EQ. 1) THEN
       H=GAUSIN(HF*(1-Q))**2
      ELSEIF(N .EQ. 2) THEN
       H=-2*LOG(1-Q)
      ELSE
       F1=Z1/N
       T=GAUSIN(Q)
       F2=SQRT(F1)*T
       IF(N .LT. 2+INT(4*ABS(T))) THEN
        S1=C1(1)
        S2=C2(1)
        S3=C3(1)
        DO 1 I = 2,7
        S1=C1(I)+S1*F2
        S2=C2(I)+S2*F2
    1   S3=C3(I)+S3*F2
        H=(S1*F1+S2)*F1+S3
       ELSE
        H=(((A(1)+A(2)*F2)*F1+(((A(3)+A(4)*F2)*F2+A(5))*F2+A(6)))*F1
     1    +(((((A(7)+A(8)*F2)*F2+A(9))*F2+A(10))*F2+A(11))*F2+A(12)))*F1
     2    +(((((A(13)*F2+A(14))*F2+A(15))*F2+A(16))*F2+A(17))*F2*F2
     3    +A(18))*F2+A(19)
       ENDIF
       H=N*H**3
      ENDIF
      CHISIN=H
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT Q = ',1P,E15.5)
  102 FORMAT('ILLEGAL DEGREE OF FREEDOM N = ',I5)
      END
