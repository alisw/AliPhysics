#include "pdf/pilot.h"
      SUBROUTINE GSXCOR(Y,NP)
C
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +       Y(NP)
      N=1
      DO 10 IX=1,20,2
         Y(N)=    (IX)/2000.0D0
      N=N+1
   10 CONTINUE
      DO 20 IX=30,200,10
         Y(N)=    (IX)/2000.0D0
           N=N+1
   20 CONTINUE
      DO 30 IX=240,1600,40
         Y(N)=    (IX)/2000.0D0
      N=N+1
   30 CONTINUE
      DO 40 IX=1625,1980,25
         Y(N)=    (IX)/2000.0D0
      N=N+1
   40 CONTINUE
      RETURN
      END
