*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION THETA1(U,Q)
C-----COMPUTATION OF THETA-SERIES. U=ARGUMENT, Q=JACOBI'S NOME.
      THETA1=2.*SIN(U)
      GO TO 1
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY THETA2
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY THETA2(U,Q)
#endif
      THETA1=2.*COS(U)
    1 IF(THETA1.EQ.0.) RETURN
      C=THETA1**2-2.
      B=THETA1*SQRT(SQRT(Q))
      A=C*B-B
      QQ=Q**2
      THETA1=QQ*A+B
      QM=QQ
      QN=QQ
      GO TO 3
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY THETA3
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY THETA3(U,Q)
#endif
      A= 2.*COS(U+U)
      GO TO 2
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY THETA4
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY THETA4(U,Q)
#endif
      A=-2.*COS(U+U)
    2 C=A
      THETA1=Q*A+1.
      QM=Q
      QN=Q
      QQ=Q**2
      B=2.
    3 DO 4 I=1,30
      QM=QQ*QM
      QN=QM*QN
      IF(QN.LT.3E-15) GO TO 5
      B=C*A-B
      THETA1=QN*B+THETA1
      QM=QQ*QM
      QN=QM*QN
      IF(QN.LT.3E-15) GO TO 5
      A=C*B-A
    4 THETA1=QN*A+THETA1
    5 RETURN
      END
