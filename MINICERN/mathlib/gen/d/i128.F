*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE I128(C,D,IX,F,G,IZ)
      REAL C(1),D(1),F(1),G(1),A(128),B(128)
      COMMON /TAB/E(253)
      IND = 0
      GOTO 20
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY D128
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY D128(C,D,IX,F,G,IZ)
#endif
      IND = 1
   20 IT = 1+IX
      JT = IX+IX
      CALL I64(C,D,JT,A,A(2),2)
      CALL I64(C(IT),D(IT),JT,B,B(2),2)
      IF(IND.EQ.1) GOTO 30
      IS = IZ
      MZ = 1
      LZ = 32*IZ
      KZ = 1+LZ
      JZ = KZ+LZ
      LZ = JZ+LZ
      GOTO 40
   30 IS = -IZ
      MZ = 32*IZ
      LZ = 1+MZ
      JZ = MZ+LZ
      KZ = JZ+MZ
      MZ = KZ+MZ
   40 F(1) = A(1)+B(1)
      G(1) = A(2)+B(2)
      F(JZ) = A(1)-B(1)
      G(JZ) = A(2)-B(2)
      F(KZ) = A(65)+B(66)
      G(KZ) = A(66)-B(65)
      F(LZ) = A(65)-B(66)
      G(LZ) = A(66)+B(65)
      DO 10 I=1,61,2
      K = I+66
      MZ = MZ+IS
      JZ = JZ+IS
      KZ = KZ+IS
      LZ = LZ+IS
      E1 = E(2*I+2)
      E2 = E(2*I+1)
      B1R = B(I+2)*E1+B(I+3)*E2
      B1I = B(I+3)*E1-B(I+2)*E2
      B2R = B(K+1)*E1-B(K)*E2
      B2I = -B(K)*E1-B(K+1)*E2
      F(MZ) = A(I+2)+B1R
      G(MZ) = A(I+3)+B1I
      F(JZ) = A(I+2)-B1R
      G(JZ) = A(I+3)-B1I
      F(KZ) = A(K)+B2R
      G(KZ) = A(K+1)+B2I
      F(LZ) = A(K)-B2R
   10 G(LZ) = A(K+1)-B2I
      RETURN
      END
