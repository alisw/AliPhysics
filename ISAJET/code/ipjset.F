#include "isajet/pilot.h"
      SUBROUTINE IPJSET
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       Initialize PJSET starting from PJETS
C-
C-   Created  14-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/primar.inc"
#include "isajet/pjets.inc"
#include "isajet/jetset.inc"
      INTEGER I,K
C----------------------------------------------------------------------
      DO 110 I=1,NJET
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*I
        JTYPE(NJSET)=IDJETS(I)
        JDCAY(NJSET)=0
        DO 115 K=1,5
115     PJSET(K,NJSET)=PJETS(K,I)
        IFRAME(I)=1
110   CONTINUE
  999 RETURN
      END
