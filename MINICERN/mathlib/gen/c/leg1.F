*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LEG1
#include "legbl.inc"
     1,SIGMA,FAC
     2,CLOGOK,CDIGAM
      DATA GAMMA/0.5772156649/
      A=-VV
      B=VV+1.0
      C=(1.0,0.0)
      CALL HYPGM(A,B,C,Z1,PP,ACC,NCVG)
      SIGMA=(0.0,0.0)
      FAC=(1.0,0.0)
      A=0.5*CLOGOK((ZZ+1.0)/(ZZ-1.0),NFRIG,1)-GAMMA-CDIGAM(B)
      QQ=A
      DO 17 L=1,50
      EL=L
      SIGMA=SIGMA+1.0/EL
      FAC=-FAC*(VV+EL)*(VV-EL+1.0)*Z1/(EL*EL)
      B=(A+SIGMA)*FAC
      QQ=QQ+B
      IF(CRIT(QQ,B,ACC))80,17,17
   17 CONTINUE
      NCVG=NCVG+2
   80 RETURN
      END
