*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION VAVDEN(RLAM)

      COMMON /GXXXC1/ WCM(0:200),AC(0:13),HC(0:8),ITYPE,NPT

      DIMENSION H(9)

      V=0
      IF(RLAM .LT. AC(0) .OR. RLAM .GT. AC(8)) GO TO 9

      IF(ITYPE .EQ. 1) THEN
       FN=1
       X=(RLAM+HC(0))*HC(1)
       H(1)=X
       H(2)=X**2-1
       DO 1 K = 2,8
       FN=FN+1
    1  H(K+1)=X*H(K)-FN*H(K-1)
       S=1+HC(7)*H(9)
       DO 2 K = 2,6
    2  S=S+HC(K)*H(K+1)
       V=HC(8)*EXP(-0.5*X**2)*MAX(S,0.)
      ELSEIF(ITYPE .EQ. 2) THEN
       X=RLAM**2
       V=AC(1)*EXP(-AC(2)*(RLAM+AC(5)*X)-
     1   AC(3)*EXP(-AC(4)*(RLAM+AC(6)*X)))
      ELSEIF(ITYPE .EQ. 3) THEN
       IF(RLAM .LT. AC(7)) THEN
        X=RLAM**2
        V=AC(1)*EXP(-AC(2)*(RLAM+AC(5)*X)-
     1         AC(3)*EXP(-AC(4)*(RLAM+AC(6)*X)))
       ELSE
        X=1/RLAM
        V=(AC(11)*X+AC(12))*X
       ENDIF
      ELSEIF(ITYPE .EQ. 4) THEN
       V=AC(13)*DENLAN(RLAM)
      ENDIF
    9 VAVDEN=V
      RETURN

      ENTRY VAVDIS(RLAM)

      IF(RLAM .LE. AC(0)) THEN
       V=0
      ELSEIF(RLAM .GE. AC(8)) THEN
       V=1
      ELSE
       X=RLAM-AC(0)
       K=X*AC(10)
       V=MIN(WCM(K)+(X-K*AC(9))*(WCM(K+1)-WCM(K))*AC(10),1.)
      ENDIF
      VAVDIS=V
      RETURN

      ENTRY VAVRND(RAN)

      K=MAX(MIN(ABS(LOCATF(WCM(0),NPT+1,RAN)),NPT)-1,0)
      X=WCM(K+1)-WCM(K)
      V=AC(0)+K*AC(9)
      IF(X .GT. 0) V=V+(RAN-WCM(K))*AC(9)/X
      VAVRND=V
      RETURN
      END
