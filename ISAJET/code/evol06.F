#include "isajet/pilot.h"
      SUBROUTINE EVOL06
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 6 (WPAIR)
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
#include "isajet/keys.inc"
#include "isajet/frame.inc"
      REAL    OFF,BP
      INTEGER I,K,J,NJSAVE,NJFINL,JTRUE
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
C
C          Add extra momenta for WPAIR
      N0PAIR=NJSET+1
      DO 130 J=1,NPAIR
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*JPAIR(J)
        JTYPE(NJSET)=IDPAIR(J)
        JDCAY(NJSET)=0
        DO 135 K=1,5
135     PJSET(K,NJSET)=PPAIR(K,J)
130   CONTINUE
      DO 140 J=1,NPAIR,2
        JET=JPAIR(J)
        JTRUE=N0PAIR+J-1
        JDCAY(N0JETS+JET-1)=JTRUE*JPACK+JTRUE+1
140   CONTINUE
      NJSAVE=NJSET
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(1.0)
C
      IF(NJSET.LT.0) RETURN
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      DO 200 I=3,NJSAVE,2
        JMATCH(I)=I+1
200   JMATCH(I+1)=I
      DO 230 I=1,2
        DO 231 K=1,5
231     FRAME(K,I)=PJSET(K,N0JETS+I-1)
        IFRAME(I)=I
230   CONTINUE
C
C          Set up and generate final state QCD parton shower.
C          Boost PJSET with -FRAME.
C
      DO 240 J=1,NJSAVE
        JET=JORIG(J)/JPACK
        IF(JET.EQ.0) JET=3
        IF(JET.GT.10) GO TO 240
        IF(IDJETS(JET).EQ.10) GO TO 240
C          Do this boost in double precision for 32-bit machines
        CALL DBOOST(-1,FRAME(1,JET),PJSET(1,J))
240   CONTINUE
C
C          Set maximum off-shell masses and JDCAY flags.
C
      NJFINL=N0PAIR
      DO 330 J=1,NPAIR
        IF(IABS(JTYPE(N0PAIR+J-1)).LT.10) THEN
          PJSET(5,N0PAIR+J-1)=PJETS(5,JPAIR(J))
          JDCAY(N0PAIR+J-1)=-1
        ENDIF
330   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END
