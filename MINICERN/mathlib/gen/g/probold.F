*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_OBSOLETE)
C     This corresponds to PROB,IF=DOUBLE and PROB64,IF=-DOUBLE
      FUNCTION PROB(CHI2,N)

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'PROB')
      PARAMETER (R1 = 1, HF = R1/2, TH = R1/3, F1 = 2*R1/9)
      PARAMETER (C1 = 1.12837 91670 95513D0)
C      PARAMETER (UP = 340)
      PARAMETER (UP = 170)

      X=HF*CHI2
      IF(N .LE. 0) THEN
       H=0
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAME,'G100.1',ERRTXT)
      ELSEIF(CHI2 .LT. 0) THEN
       H=0
       WRITE(ERRTXT,102) CHI2
       CALL MTLPRT(NAME,'G100.2',ERRTXT)
      ELSEIF(CHI2 .GT. UP) THEN
       H=0
      ELSEIF(N .GT. 100) THEN
       S=R1/N
       T=F1*S
       H=HF*ERFC(((CHI2*S)**TH-(1-T))/SQRT(2*T))
      ELSEIF(N .EQ. 1) THEN
       H=ERFC(SQRT(X))
      ELSE
       S=1
       T=1
       M=N/2
       IF(2*M .EQ. N) THEN
        DO 1 I = 1,M-1
        T=X*T/I
    1   S=S+T
        H=S*EXP(-X)
       ELSE
        DO 2 I=1,M-1
        T=T*CHI2/(2*I+1)
    2   S=S+T
        W=SQRT(X)
        H=C1*W*S*EXP(-X)+ERFC(W)
       ENDIF
      ENDIF
      PROB=H
      RETURN
  101 FORMAT('N = ',I6,' < 1')
  102 FORMAT('X = ',1P,E20.10,' < 0')
      END
#endif
