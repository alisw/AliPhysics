#include "isajet/pilot.h"
      SUBROUTINE GETPT(PT0,PTMEAN)
C          GENERATE PT WITH 1/(1+B*PT**2)**4 DISTRIBUTION
C          (APPROXIMATELY AN EXPONENTIAL FOR PT < 2 GEV.)
C          CON1=16/(3*PI)
C          CON2=-1/3
#include "isajet/itapes.inc"
      DATA CON1/1.697652726/,CON2/-.3333333333/
      R=RANF()
      ARG=AMAX1(R**CON2-1.,0.)
      PT0=PTMEAN*CON1*SQRT(ARG)
      RETURN
      END
