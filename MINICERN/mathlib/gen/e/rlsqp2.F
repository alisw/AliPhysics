*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RLSQP2(N,X,Y,A0,A1,A2,SD,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y,A0,A1,A2,SD
 
      DIMENSION X(*),Y(*)
 
      PARAMETER (R0 = 0)
 
      A0=0
      A1=0
      A2=0
      SD=0
      IF(N .LE. 2) THEN
       IFAIL=1
      ELSE
       FN=N
       XM=0
       DO 1 K = 1,N
       XM=XM+X(K)
    1  CONTINUE
       XM=XM/FN
       SX=0
       SXX=0
       SXXX=0
       SXXXX=0
       SY=0
       SYY=0
       SXY=0
       SXXY=0
       DO 2 K = 1,N
       XK=X(K)-XM
       YK=Y(K)
       XK2=XK**2
       SX=SX+XK
       SXX=SXX+XK2
       SXXX=SXXX+XK2*XK
       SXXXX=SXXXX+XK2**2
       SY=SY+YK
       SYY=SYY+YK**2
       SXY=SXY+XK*YK
       SXXY=SXXY+XK2*YK
    2  CONTINUE
       DET=(FN*SXXXX-SXX**2)*SXX-FN*SXXX**2
       IF(DET .GT. 0) THEN
        D2=(SXX*(FN*SXXY-SXX*SY)-FN*SXXX*SXY)/DET
        D1=(SXY-SXXX*D2)/SXX
        D0=(SY-SXX*D2)/FN
        IFAIL=0
       ELSE
        IFAIL=-1
       ENDIF
      ENDIF
      IF(IFAIL .EQ. 0 .AND. N .GT. 3)
     1   SD=SQRT(MAX(R0,SYY-D0*SY-D1*SXY-D2*SXXY)/(N-3))
      A0=D0+XM*(XM*D2-D1)
      A1=D1-2*XM*D2
      A2=D2
      RETURN
      END
