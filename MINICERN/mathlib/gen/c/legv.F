*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LEGV
#include "legbl.inc"
     1,F1,F2,SSRZ
     2,CLOGOK,CSQRTK,RGAM
      A=0.5
      C=VV+1.5
      CALL HYPGM(A,A,C,ZZ1,F1,ACC,NCVG)
      CALL HYPGM(A,A,C,ZZ2,F2,ACC,NCV )
      NCVG=NCVG+2*NCV
C*UL6 F1=F1*EXP((VV+0.5)*CLOGOK(ZZ +SRZ,NFRIG,3))
      F1=F1*EXP((VV+0.5)*CLOGOK(ZZ +SRZ,NFRIG,3))
      F2=F2*EXP((VV+0.5)*CLOGOK(ZZ -SRZ,-NFRIG,3))
      A=1.5
      B=1.0
      SSRZ=CSQRTK(2.0*SRZ,NFRIG,2)*RGAM(VV,A,B)
      SGN=1.0
      IF(AIMAG(ZZ))8,14,12
    8 SGN=-1.0
      GOTO 12
   14 SGN=SIGN(1.0,(NFRIG+0.5)*REAL(ZZ))
   12 PP=(F1+SGN*U*F2)/(PISR*SSRZ)
      IF(NFRIG)11,10,11
   10 QQ=0.5*PISR*(F2+SGN*U*F1)/SSRZ
      GOTO 80
   11 QQ=PISR*F2/SSRZ
   80 RETURN
      END
