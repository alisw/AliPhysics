*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DLSQP1(N,X,Y,A0,A1,SD,IFAIL)
#if !defined(CERNLIB_DOUBLE)
      CHARACTER*6 NAME
      NAME = 'DLSQP1'
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
        A0=(SY*SXX-SXY*SX)/DET
        A1=(FN*SXY-SX*SY)/DET
        IFAIL=0
       ELSE
        IFAIL=-1
       ENDIF
      ENDIF
      IF(IFAIL .EQ. 0 .AND. N .GT. 2)
     1   SD=SQRT(MAX(R0,SYY-A0*SY-A1*SXY)/(N-2))
      A0=A0-A1*XM
      RETURN
      END
#endif
