#include "isajet/pilot.h"
      SUBROUTINE LSTSQ(X,Y,NPT,A,B)
C
C          DO LEAST SQUARE FIT TO A STRAIGHT LINE Y=A+B*X
C
#include "isajet/itapes.inc"
      DIMENSION X(NPT),Y(NPT)
      SUM1=0
      SUM2=0
      SUM3=0
      SUM4=0
      DO 1 I=1,NPT
      SUM1=SUM1+X(I)
      SUM2=SUM2+Y(I)
      SUM3=SUM3+X(I)**2
      SUM4=SUM4+X(I)*Y(I)
    1 CONTINUE
      B=(SUM2*SUM1-SUM4*NPT)/(SUM1**2-SUM3*NPT)
      A=(SUM2-B*SUM1)/NPT
      RETURN
      END
