#include "isajet/pilot.h"
      SUBROUTINE SSNORM(ID)
C-----------------------------------------------------------------------
C          Normalize branching ratios for ID
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
C
      INTEGER ID,I
      REAL GAMSUM
C
      GAMSUM=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).EQ.ID) GAMSUM=GAMSUM+GSSMOD(I)
100   CONTINUE
      IF(GAMSUM.EQ.0) RETURN
      DO 200 I=1,NSSMOD
        IF(ISSMOD(I).EQ.ID) BSSMOD(I)=GSSMOD(I)/GAMSUM
200   CONTINUE
      RETURN
      END 
