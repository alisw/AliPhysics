*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:24  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DLSQPM(N,X,Y,M,A,SD,IFAIL)
#if !defined(CERNLIB_DOUBLE)
      CHARACTER*6 NAME
      NAME = 'DLSQPM'
      CALL MTLPRT(NAME,'E201',
     +'not available on this machine - see documentation')
      RETURN
      END
#endif
#if defined(CERNLIB_DOUBLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      PARAMETER (IDIM = 21, R0 = 0)
 
      DIMENSION  X(*),Y(*),A(0:*),B(IDIM,IDIM),XY(0:IDIM)
 
      M1=M+1
      IFAIL=0
      IF(N .LE. 1 .OR. M .LT. 0 .OR. M1 .GT. IDIM .OR. M1 .GT. N) THEN
       IFAIL=1
      ELSEIF(M .EQ. 0) THEN
       XY(0)=DVSUM(N,Y(1),Y(2))
       A(0)=XY(0)/N
       SYY=DVMPY(N,Y(1),Y(2),Y(1),Y(2))
      ELSE
       DO 11 J = 1,M1
       A(J-1)=0
       B(J,1)=0
       B(M1,J)=0
   11  CONTINUE
       B(1,1)=N
       SYY=0
       DO 4 K = 1,N
       XK=X(K)
       YK=Y(K)
       SYY=SYY+Y(K)**2
       P=1
       A(0)=A(0)+YK
       DO 2 J = 2,M1
       P=P*XK
       B(J,1)=B(J,1)+P
       A(J-1)=A(J-1)+P*YK
    2  CONTINUE
       DO 3 J = 2,M1
       P=P*XK
       B(M1,J)=B(M1,J)+P
    3  CONTINUE
    4  CONTINUE
       DO 5 I = 2,M
       DO 5 K = I,M1
       B(K-1,I)=B(K,I-1)
    5  CONTINUE
       DO 6 J = 0,M
       XY(J)=A(J)
    6  CONTINUE
       CALL DSEQN(M1,B,IDIM,IFAIL,1,A)
      ENDIF
      SD=0
      IF(IFAIL .EQ. 0) THEN
       IF(N .GT. M1) THEN
        SD=SYY
        DO 7 J = 0,M
        SD=SD-A(J)*XY(J)
    7   CONTINUE
        SD=SQRT(MAX(R0,SD)/(N-M1))
       ENDIF
      ELSE
       CALL DVSET(M1,R0,A(0),A(1))
       M=0
      ENDIF
      RETURN
      END
#endif
