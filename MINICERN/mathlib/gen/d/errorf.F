*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:21  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE ERRORF (M)
C-----FINAL ERROR ANALYSIS FOR CHISQ.
#include "d510si.inc"
#include "d510ui.inc"
#include "d510uo.inc"
      EQUIVALENCE(NED(1),NS),(NED(2),NP)
      WRITE(6,3)
      IF (NP.GT.7) WRITE(6,4)
      DO 2 J=1,NS
      K1=NP*(J-1)+1
      FUN=EXDA(K1)
      SIG=EXDA(K1+1)
      NX=NP-2
      DO 1 K=1,NX
      K2=K1+1+K
      X(K)=EXDA(K2)
 1    CONTINUE
C-----CALCULATE DY/DA"S FOR SCAL
      CALL ARITHM (Y)
      CALL SCAL (M,ER)
      XI=(FUN-Y)/SIG
      XI=0.5*XI*XI
      ER=SQRT(ER)
      ERROR(J)=ER
      WRITE(6,5) J,Y,ER,XI,(X(I),I=1,NX)
 2    CONTINUE
      RETURN
C
 3    FORMAT('1','ERROR ANALYSIS FOR EACH DATA POINT USING FINAL PARA',
     +'METER VALUES (CHISQ ONLY)'////11X,'POINT',8X,'FITTED',9X,
     +'Y STANDARD',6X,'CONTRIBUTION',4X,'X CO-ORDINATES'/12X,'NO.',8X,
     +'Y  VALUE',8X,'DEVIATION',7X,'TO CHISQ/2',6X,'OF DATA POINT'//)
 4    FORMAT('+',86X,'(PRINTED IN ROW ORDER)')
 5    FORMAT(11X,I3,2X,7(5X,E12.5),(/67X,4(5X,E12.5)))
      END
