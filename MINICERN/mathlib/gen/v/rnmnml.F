*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RNMNML(N,NSUM,PCUM,NVEC,IERR)
 
      PARAMETER (NRV = 100)
      DIMENSION PCUM(*),NVEC(*),RV(NRV)
 
      IERR=2
      IF(PCUM(N) .NE. 1) GO TO 9
      IERR=1
      DO 1 I = 2,N
      IF(PCUM(I) .LT. PCUM(I-1)) GO TO 9
    1 CONTINUE
      IERR=0
      DO 2 K = 1,N
    2 NVEC(K)=0
      J0=NSUM/NRV+1
      I0=NRV
      DO 3 J = 1,J0
      IF(J .EQ. J0) I0=MOD(NSUM,NRV)
      CALL RANLUX(RV,I0)
      DO 4 I = 1,I0
      DO 5 K = 1,N
      IF(PCUM(K) .GE. RV(I)) THEN
       NVEC(K)=NVEC(K)+1
       GO TO 4
      ENDIF
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    9 RETURN
      END
