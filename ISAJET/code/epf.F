#include "isajet/pilot.h"
      FUNCTION EPF(A,B,C,D)
C          CALCULATE TOTALLY ANTISYMMETRIC TENSOR EPSILON CONTRACTED
C          WITH FOUR 4-VECTORS.
#include "isajet/itapes.inc"
      DIMENSION A(4),B(4),C(4),D(4)
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION EPF
      DOUBLE PRECISION A,B,C,D,CD,BCD
#endif
      CD(I,J)=C(I)*D(J)-C(J)*D(I)
      BCD(I,J,K)=B(I)*CD(J,K)-B(J)*CD(I,K)+B(K)*CD(I,J)
      EPF=A(1)*BCD(2,3,4)-A(2)*BCD(1,3,4)+A(3)*BCD(1,2,4)
     1-A(4)*BCD(1,2,3)
      RETURN
      END
