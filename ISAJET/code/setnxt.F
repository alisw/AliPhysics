#include "isajet/pilot.h"
      SUBROUTINE SETNXT
C
C            RESET LIMITS BEFORE NEXT SET
C
#include "isajet/itapes.inc"
#include "isajet/lstprt.inc"
#include "isajet/totals.inc"
#include "isajet/dylim.inc"
#include "isajet/jetlim.inc"
#include "isajet/primar.inc"
#include "isajet/jetset.inc"
#include "isajet/partcl.inc"
      DATA UNDEF/-1.E9/
      DO 1 I=1,36
      IF(SETLMJ(I)) BLIMS(I)=UNDEF
    1 CONTINUE
      DO 2 I=1,12
      IF(SETLMQ(I)) BLIM1(I)=UNDEF
    2 CONTINUE
C          RESET /TOTALS/
      NKINPT=0
      NWGEN=0
      NKEEP=0
      SUMWT=0.
C          RESET /LSTPRT/
      LSTPRT=0
C          RESET NJSET AND NPTCL
      NJSET=0
      NPTCL=0
      NPAIR=0
      RETURN
      END
