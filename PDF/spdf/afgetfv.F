*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:31  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:59  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       FUNCTION AFGETFV(X,FVL)
C  NOUVEAU PROGRAMME D'INTERPOLATION UTILISANT UNE ROUTINE DE MATH. RECIPES
       DIMENSION FVL(32)
#include "pdf/expdp.inc"
     +                XI(32),WI(32),XX(33)
       COMMON/W5051I9/XI,WI,XX,NTERMS
       DIMENSION A(4),B(4)
       N=4
       EPS=1.E-7
       XAM=XX(1)-EPS
       XAP=XX(1)+EPS
C      IF(X.LT.XAM) PRINT*,' X = ',X
       IF(X.GT.XAM.AND.X.LT.XAP) GO TO 50
       GO TO 80
   50  Y=FVL(1)
       GO TO 77
   80  IF(X.LT.XX(2)) GO TO 51
       IF(X.GT.XX(30)) GO TO 61
       DO 1 I=3,30
       IF(X.GT.XX(I)) GO TO 1
       A(1)=XX(I-2)
       A(2)=XX(I-1)
       A(3)=XX(I)
       A(4)=XX(I+1)
       B(1)=FVL(I-2)
       B(2)=FVL(I-1)
       B(3)=FVL(I)
       B(4)=FVL(I+1)
       GO TO 70
   1   CONTINUE
  61   A(1)=XX(29)
       A(2)=XX(30)
       A(3)=XX(31)
       A(4)=XX(32)
       B(1)=FVL(29)
       B(2)=FVL(30)
       B(3)=FVL(31)
       B(4)=FVL(32)
       GO  TO 70
  51   A(1)=XX(1)
       A(2)=XX(2)
       A(3)=XX(3)
       A(4)=XX(4)
       B(1)=FVL(1)
       B(2)=FVL(2)
       B(3)=FVL(3)
       B(4)=FVL(4)
C 70   IF(X.GT..2.AND.X.LT..8) THEN
C            CALL AFPOLIN(A,B,N,X,Y,DY)
C      ELSE
C            CALL AFRATIN(A,B,N,X,Y,DY)
C      ENDIF
  70   CONTINUE
             CALL AFPOLIN(A,B,N,X,Y,DY)
  77   AFGETFV=Y
       RETURN
       END
