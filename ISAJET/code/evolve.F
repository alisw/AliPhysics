#include "isajet/pilot.h"
      SUBROUTINE EVOLVE
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Call for each process a subroutine to set up
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
#include "isajet/pinits.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/jwork2.inc"
#include "isajet/keys.inc"
#include "isajet/frame.inc"
      REAL BP,PINCOM
      INTEGER I,K,J,JJET,IFR
C----------------------------------------------------------------------
C          Initialize
      NJSET=0
      N0JETS=0
      N0W=0
      N0PAIR=0
C
C          Copy momenta from /PINITS/ to /JETSET/
      IF(.NOT.KEYS(2)) THEN
        DO 100 I=1,2
          NJSET=NJSET+1
          JORIG(NJSET)=JPACK*(10+I)
          JTYPE(NJSET)=IDINIT(I)
          JDCAY(NJSET)=JPACK*I+I
          DO 105 K=1,5
105       PJSET(K,NJSET)=PINITS(K,I)
100     CONTINUE
      ENDIF
C
C       Handle each process separately
C
      IF(KEYS(1).OR.KEYS(8)) THEN
        CALL EVOL01
      ELSEIF(KEYS(2)) THEN
        CALL EVOL02
      ELSEIF(KEYS(3)) THEN
        CALL EVOL03
      ELSEIF(KEYS(5)) THEN
        CALL EVOL05
      ELSEIF(KEYS(6).OR.KEYS(10)) THEN
        CALL EVOL06
      ELSEIF(KEYS(7).OR.KEYS(9)) THEN
        CALL EVOL07
      ELSEIF(KEYS(11)) THEN
        CALL EVOL11
      ELSEIF(KEYS(12)) THEN
        CALL EVOL01
      ENDIF
C
      IF(NJSET.LT.0) RETURN
C
C          Boost /JETSET/ partons back to PP COM
C
      DO 500 J=1,NJSET
        JJET=JORIG(J)/JPACK
        IF ( JJET.EQ.0 ) THEN
          IFR=1
        ELSE
          IF(JJET.GT.10) GO TO 500
          IF(IDJETS(JJET).EQ.10.AND.KEYS(6)) GO TO 500
          IFR=IFRAME(JJET)
        ENDIF
        BP=0.
        DO 505 K=1,3
505     BP=BP+FRAME(K,IFR)*PJSET(K,J)
        BP=BP/FRAME(5,IFR)
        DO 510 K=1,3
510     PJSET(K,J)=PJSET(K,J)+FRAME(K,IFR)*PJSET(4,J)/FRAME(5,IFR)
     1  +FRAME(K,IFR)*BP/(FRAME(4,IFR)+FRAME(5,IFR))
        PJSET(4,J)=FRAME(4,IFR)*PJSET(4,J)/FRAME(5,IFR)+BP
500   CONTINUE
C
C          Reset PBEAM
      DO 530 J=1,NJSET
        IF(JDCAY(J).EQ.JPACK*J+J) THEN
          JJET=JORIG(J)/JPACK-10
          PINCOM=.5*(PJSET(4,J)+ABS(PJSET(3,J)))
          PBEAM(JJET)=HALFE-PINCOM
        ENDIF
530   CONTINUE
C
C          Check for zero energy partons
      CALL IRMOV0
C
      RETURN
      END
