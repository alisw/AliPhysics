#include "isajet/pilot.h"
      SUBROUTINE EVOL03
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 3 (DRELLYAN)
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C-   Created  13-AUG-1991   Frank E. Paige,Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/pjets.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/jwork2.inc"
#include "isajet/q1q2.inc"
#include "isajet/frame.inc"
#include "isajet/wcon.inc"
      REAL    EVOLMS,BP
      INTEGER I,K,J,NJFINL
C----------------------------------------------------------------------
C
C          Add W momentum and recoil jets
      N0JETS=NJSET+1
      IF(.NOT.STDDY) THEN
        DO 101 I=3,NJET
          NJSET=NJSET+1
          JORIG(NJSET)=JPACK*I
          JTYPE(NJSET)=IDJETS(I)
          JDCAY(NJSET)=0
          DO 105 K=1,5
105       PJSET(K,NJSET)=PJETS(K,I)
          IFRAME(I)=1
101     CONTINUE
        NJSET=NJSET+1
        N0W=NJSET
        JORIG(NJSET)=0
        JTYPE(NJSET)=IDENTW
        JDCAY(NJSET)=(N0W+1)*JPACK+N0W+2
        DO 120 K=1,5
120     PJSET(K,NJSET)=QWJET(K)
      ENDIF
C
C          Add W decays
      DO 110 I=1,2
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*I
        JTYPE(NJSET)=IDJETS(I)
        JDCAY(NJSET)=0
        DO 115 K=1,5
115     PJSET(K,NJSET)=PJETS(K,I)
        IFRAME(I)=2
        IF(STDDY) IFRAME(I)=1
110   CONTINUE
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(WFUDGE)
C
      IF(NJSET.LT.0) RETURN
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      IF(STDDY) THEN
        CALL IFRAMS(3,4,1,.FALSE.)
      ELSE
        CALL IFRAMS(N0W+1,N0W+2,2,.FALSE.)
        CALL IFRAMS(N0JETS,N0W,1,.FALSE.)          
      ENDIF
C
C          Set maximum off-shell masses and JDCAY flags.
C
      IF(STDDY) THEN
        NJFINL=3
        DO 310 J=3,4
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=QMW
            JDCAY(J)=-1
          ENDIF
310     CONTINUE
      ELSE
        NJFINL=N0JETS
        DO 320 J=N0W+1,N0W+2
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=QMW
            JDCAY(J)=-1
          ENDIF
320     CONTINUE
C          Need fudge factor for DRELLYAN
        DO 321 J=N0JETS,N0W
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=EVOLMS(J,WFUDGE)
            JDCAY(J)=-1
          ENDIF
321     CONTINUE
      ENDIF
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
C          Reset FRAME using W momentum modified by evolution
      IF(.NOT.STDDY) THEN
        BP=0.
        DO 400 K=1,3
400     BP=BP+FRAME(K,1)*PJSET(K,N0W)
        BP=BP/FRAME(5,1)
        DO 410 K=1,3
          FRAME(K,2)=PJSET(K,N0W)+FRAME(K,1)*PJSET(4,N0W)/FRAME(5,1)
     $    +FRAME(K,1)*BP/(FRAME(4,1)+FRAME(5,1))
410     CONTINUE
        FRAME(4,2)=FRAME(4,1)*PJSET(4,N0W)/FRAME(5,1)+BP
      ENDIF
C
      RETURN
      END
