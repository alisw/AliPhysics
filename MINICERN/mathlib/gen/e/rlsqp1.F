*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RLSQP1(N,X,Y,A0,A1,SD,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y,A0,A1,SD
 
      DIMENSION X(*),Y(*)

      PARAMETER (R0 = 0)

      A0=0
      A1=0
      SD=0
      IF(N .LE. 1) THEN
       IFAIL=1
      ELSE
       FN=N
       XM=0
       DO 1 K = 1,N
       XM=XM+X(K)
    1  CONTINUE
       XM=XM/FN
       SX=0
       SY=0
       SXX=0
       SYY=0
       SXY=0
       DO 2 K = 1,N
       XK=X(K)-XM
       YK=Y(K)
       SX=SX+XK
       SY=SY+YK
       SXX=SXX+XK**2
       SYY=SYY+YK**2
       SXY=SXY+XK*YK
    2  CONTINUE
       DET=FN*SXX-SX**2
       IF(DET .GT. 0) THEN
        D0=(SY*SXX-SXY*SX)/DET
        D1=(FN*SXY-SX*SY)/DET
        IFAIL=0
       ELSE
        IFAIL=-1
       ENDIF
      ENDIF
      IF(IFAIL .EQ. 0) THEN
       IF(N .GT. 2) SD=SQRT(MAX(R0,SYY-D0*SY-D1*SXY)/(N-2))
       A0=D0-D1*XM
       A1=D1
      ENDIF
      RETURN
      END
