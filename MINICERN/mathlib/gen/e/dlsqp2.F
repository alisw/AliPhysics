*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DLSQP2(N,X,Y,A0,A1,A2,SD,IFAIL)
#if !defined(CERNLIB_DOUBLE)
      CHARACTER*6 NAME
      NAME = 'DLSQP2'
      CALL MTLPRT(NAME,'E201',
     +'not available on this machine - see documentation')
      RETURN
      END
#endif
#if defined(CERNLIB_DOUBLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
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
        A2=(SXX*(FN*SXXY-SXX*SY)-FN*SXXX*SXY)/DET
        A1=(SXY-SXXX*A2)/SXX
        A0=(SY-SXX*A2)/FN
        IFAIL=0
       ELSE
        IFAIL=-1
       ENDIF
      ENDIF
      IF(IFAIL .EQ. 0 .AND. N .GT. 3)
     1   SD=SQRT(MAX(R0,SYY-A0*SY-A1*SXY-A2*SXXY)/(N-3))
      A0=A0+XM*(XM*A2-A1)
      A1=A1-2*XM*A2
      RETURN
      END
#endif
