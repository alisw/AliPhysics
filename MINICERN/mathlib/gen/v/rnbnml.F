*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RNBNML(M,P,N,IERR)
 
      PARAMETER (NRV = 100)
      DIMENSION RV(NRV)
 
      IERR=1
      IF(0 .LT. P .AND. P .LT. 1) THEN
       IERR=0
       N=0
       J0=M/NRV+1
       I0=NRV
       DO 1 J = 1,J0
       IF(J .EQ. J0) I0=MOD(M,NRV)
       CALL RANLUX(RV,I0)
       DO 1 I = 1,I0
       IF(RV(I) .LE. P) N=N+1
    1  CONTINUE
      ENDIF
      RETURN
      END
