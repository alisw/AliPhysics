*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION RNGAMA(P)
      EXTERNAL RANLUX
 
      DIMENSION STOR(15)
 
      H=0
      IF(P .GT. 15) THEN
    1  CALL RNORMX(A,1,RANLUX)
       H=P*(1-1/(9*P)+A/(3*SQRT(P)))**3
       IF(H .LE. 0) GO TO 1
      ELSE
       M=P
       F=P-M
       IF(M .GT. 0) THEN
        X=1
        CALL RANLUX(STOR,M)
        DO 2 I = 1,M
    2   X=X*STOR(I)
        H=-LOG(X)
       ENDIF
       IF(F .GE. 0.00001) THEN
        CALL RANLUX(X,1)
        X1=-LOG(X)
        IF(F .GE. 0.9999) THEN
         H=H+X1
        ELSE
         CALL RANLUX(X,1)
    3    WLOG=LOG(X)/F
         IF(WLOG .GE. -100) THEN
          W1=EXP(WLOG)
          CALL RANLUX(X,1)
          WLOG=LOG(X)/(1-F)
          IF(WLOG .LT. -100) THEN
           H=H+X1
          ELSE
           W=W1+EXP(WLOG)
           IF(W .GT. 1) GO TO 3
           H=H+X1*W1/W
          ENDIF
         ENDIF
        ENDIF
       ENDIF
      ENDIF
      RNGAMA=H
      RETURN
      END
