*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FFTRC(C,D,IX,F,G,IZ)
C
C         D701  FFTRC  FAST FOURIER TRANSFORM, REAL OR COMPLEX
C         B. FORNBERG  NOVEMBER 1973
C
      COMPLEX A(16),T
      DIMENSION E(32),C(1),D(1),F(1),G(1),S(2),NC(7)
      COMMON /TAB/U(126),V(127)
      EQUIVALENCE (A(1),E(1)),(S(1),T),(S(1),S1),(S(2),S2),(U(31),C3),(U
     *(63),C7),(U(32),C9)
      DATA NC/17,9,25,5,21,13,29/
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY I16
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY I16(C,D,IX,F,G,IZ)
#endif
      IS = IX
      I = 1
      GOTO 40
C
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
      ENTRY D16
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
      ENTRY D16(C,D,IX,F,G,IZ)
#endif
      IS = -IX
      I = 1+16*IX
   40 J = 1+8*IX
      CALL D701BD
      E(3) = C(1)-C(J)
      E(1) = C(1)+C(J)
      E(4) = D(1)-D(J)
      E(2) = D(1)+D(J)
      DO 10 L=1,7
      I = I+IS
      J = J+IS
      N = NC(L)
      E(N+2) = C(I)-C(J)
      E(N) = C(I)+C(J)
      E(N+3) = D(I)-D(J)
   10 E(N+1) = D(I)+D(J)
      DO 20 L=1,13,4
      T = A(L+2)
      A(L+2) = A(L)-T
   20 A(L) = A(L)+T
      DO 30 L=2,14,4
      S1 = E(2*L+4)
      S2 = -E(2*L+3)
      A(L+2) = A(L)-T
   30 A(L) = A(L)+T
      T = A(5)
      A(5) = A(1)-T
      A(1) = A(1)+T
      T = A(13)
      A(13) = A(9)-T
      A(9) = A(9)+T
      S2 = -C7*(E(11)-E(12))
      S1 = C7*(E(11)+E(12))
      A(6) = A(2)-T
      A(2) = A(2)+T
      S2 = -C7*(E(27)-E(28))
      S1 = C7*(E(27)+E(28))
      A(14) = A(10)-T
      A(10) = A(10)+T
      S1 = E(14)
      S2 = -E(13)
      A(7) = A(3)-T
      A(3) = A(3)+T
      S1 = E(30)
      S2 = -E(29)
      A(15) = A(11)-T
      A(11) = A(11)+T
      S2 = -C7*(E(15)+E(16))
      S1 = -C7*(E(15)-E(16))
      A(8) = A(4)-T
      A(4) = A(4)+T
      S1 = -C7*(E(31)-E(32))
      S2 = -C7*(E(31)+E(32))
      A(16) = A(12)-T
      A(12) = A(12)+T
      I = 1
      J = 1+8*IZ
      F(J) = E(1)-E(17)
      G(J) = E(2)-E(18)
      F(I) = E(1)+E(17)
      G(I) = E(2)+E(18)
      I = I+IZ
      J = J+IZ
      S1 = C9*E(19)+C3*E(20)
      S2 = C9*E(20)-C3*E(19)
      F(J) = E( 3)-S1
      G(J) = E( 4)-S2
      F(I) = E( 3)+S1
      G(I) = E( 4)+S2
      I = I+IZ
      J = J+IZ
      S1 = C7*(E(21)+E(22))
      S2 = -C7*(E(21)-E(22))
      F(J) = E( 5)-S1
      G(J) = E( 6)-S2
      F(I) = E( 5)+S1
      G(I) = E( 6)+S2
      I = I+IZ
      J = J+IZ
      S1 = C3*E(23)+C9*E(24)
      S2 = C3*E(24)-C9*E(23)
      F(J) = E( 7)-S1
      G(J) = E( 8)-S2
      F(I) = E( 7)+S1
      G(I) = E( 8)+S2
      I = I+IZ
      J = J+IZ
      F(J) = E(9)-E(26)
      G(J) = E(10)+E(25)
      F(I) = E(9)+E(26)
      G(I) = E(10)-E(25)
      I = I+IZ
      J = J+IZ
      S1 = -C3*E(27)+C9*E(28)
      S2 = -C3*E(28)-C9*E(27)
      F(J) = E(11)-S1
      G(J) = E(12)-S2
      F(I) = E(11)+S1
      G(I) = E(12)+S2
      I = I+IZ
      J = J+IZ
      S1 = -C7*(E(29)-E(30))
      S2 = -C7*(E(29)+E(30))
      F(J) = E(13)-S1
      G(J) = E(14)-S2
      F(I) = E(13)+S1
      G(I) = E(14)+S2
      I = I+IZ
      J = J+IZ
      S1 = -C9*E(31)+C3*E(32)
      S2 = -C9*E(32)-C3*E(31)
      F(J) = E(15)-S1
      G(J) = E(16)-S2
      F(I) = E(15)+S1
      G(I) = E(16)+S2
      RETURN
      END
