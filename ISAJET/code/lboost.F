#include "isajet/pilot.h"
      SUBROUTINE LBOOST(PREST,N,P1,P2)
C
C            BOOST 4-VECTORS P1 TO PREST REST FRAME
C            PUT RESULTING 4-VECTORS IN P2
C
#include "isajet/itapes.inc"
      DIMENSION PREST(4),P1(4,N),P2(4,N)
      DO 1 I=1,N
      WCN=SQRT(PREST(4)**2-PREST(1)**2-PREST(2)**2-PREST(3)**2)
      II=(I-1)*4
      P2(4,I)=(P1(4,I)*PREST(4)-P1(1,I)*PREST(1)-P1(2,I)*PREST(2)
     1         -P1(3,I)*PREST(3))/WCN
      FACT=(P2(4,I)+P1(4,I))/(WCN+PREST(4))
      DO 2 K=1,3
    2 P2(K,I)=P1(K,I)-FACT*PREST(K)
    1 CONTINUE
      RETURN
      END
