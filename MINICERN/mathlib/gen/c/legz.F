*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LEGZ
#include "legbl.inc"
     1,F1,F2,ZV
     2,CLOGOK,RGAM
      GR=CLOGOK(2.0*ZZ,NFRIG,2)
      A=1.5
      B=1.0
      ZV=RGAM(VV,A,B)*EXP(VV*GR)
      A=VV/2.0+1.0
      B=VV/2.0+0.5
      C=VV+1.5
      ACCC=ACC/100.
      CALL HYPGM(A,B,C,Z2,F1,ACCC,NCVG)
      F1=F1/(2.0*ZZ*ZV)
      QQ=PISR*F1
      IF(ABS(CVV).GT.0.001) GOTO 9
C   TRDZ EXPECTS GR=CLOGOK(2.0*ZZ,NFRIG,2) BUT DESTROYS CONTENTS
      CALL TRDZ
      GOTO 80
    9 A=-VV/2.0
      B=(1.0-VV)/2.0
      C=0.5-VV
      CALL HYPGM(A,B,C,Z2,F2,ACCC,NCV )
      NCVG=NCVG+2*NCV
      F2=F2*ZV/(VV+0.5)
      PP=(F1*SVV/CVV+F2)/PISR
   80 RETURN
      END
