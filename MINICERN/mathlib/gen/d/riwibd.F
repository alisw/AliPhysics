*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:16  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE RIWIBD
C   BLOCK DATA SUBSTITUTE
#include "gen/imp64.inc"
      COMMON/STORE/XA(11),XB(11),XC(11),XD(11),MA(11),MB(11),MC(11)
      COMMON/STORE1/R(10000),LR
      COMMON/OPTION/IPRINT,ICONV,IRESET
      COMMON/RANDOM/NSHOTS
      COMMON/INTERN/FACTOR,ALFA,BETA,GAMMA,DELTA,LEVEL,NMIN
      DATA INIT/0/
      IF(INIT.EQ.1) RETURN
      INIT=1
      MA(1)=0
      LR=10000
      IPRINT=1
      ICONV=1
      IRESET=0
      NSHOTS=2
      FACTOR=1.65
      LEVEL=90
      ALFA=0.3
      BETA=0.3
      GAMMA=0.3
      DELTA=.7
      NMIN=2
      RETURN
      END
