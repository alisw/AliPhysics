#include "isajet/pilot.h"
      SUBROUTINE HIGGS
C
C          FINISH HIGGS GENERATION STARTED BY DRLLYN FOR DECAY
C          HIGGS --> W W.
C
C          VER 7.14: TEST BOTH JET1 AND JET2 FOR W,Z FOR SAFETY
C
#include "isajet/itapes.inc"
#include "isajet/qcdpar.inc"
#include "isajet/jetpar.inc"
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/jetsig.inc"
#include "isajet/qsave.inc"
#include "isajet/wcon.inc"
#include "isajet/const.inc"
#include "isajet/hcon.inc"
C
      DIMENSION X(2)
      EQUIVALENCE (X(1),X1)
C
      IDABS1=IABS(IDJETS(1))
      IDABS2=IABS(IDJETS(2))
      IF(IDABS1.NE.80.AND.IDABS1.NE.90.AND.
     $IDABS2.NE.80.AND.IDABS2.NE.90) THEN
        NPAIR=0
        DO 100 I=1,4
        IDPAIR(I)=0
        JPAIR(I)=0
        DO 110 K=1,5
110     PPAIR(K,I)=0.
100     CONTINUE
      ELSE
        CALL WPAIR
      ENDIF
      RETURN
      END
