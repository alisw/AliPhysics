#include "isajet/pilot.h"
      SUBROUTINE EVOL07
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 7 (HIGGS)
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
#include "isajet/frame.inc"
      REAL    EVOLMS,BP
      INTEGER I,K,J,NJSAVE,NJFINL,JTRUE
      DOUBLE PRECISION DPASS(5),DSUM(5)
      INTEGER IDABS1,IDABS2
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
C
C          Add extra momenta for WPAIR
      IDABS1=IABS(IDJETS(1))
      IDABS2=IABS(IDJETS(2))
      IF(IDABS1.EQ.80.OR.IDABS1.EQ.90.OR.IDABS2.EQ.80.OR.
     $IDABS2.EQ.90) THEN
        N0PAIR=NJSET+1
        DO 130 J=1,NPAIR
          NJSET=NJSET+1
          JORIG(NJSET)=JPACK*JPAIR(J)
          JTYPE(NJSET)=IDPAIR(J)
          JDCAY(NJSET)=0
          DO 135 K=1,5
135       PJSET(K,NJSET)=PPAIR(K,J)
130     CONTINUE
        DO 140 J=1,NPAIR,2
          JET=JPAIR(J)
          JTRUE=N0PAIR+J-1
          JDCAY(N0JETS+JET-1)=JTRUE*JPACK+JTRUE+1
140     CONTINUE
      ENDIF
      NJSAVE=NJSET
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      IF(IABS(IDINIT(1)).LT.80) THEN
        CALL ISTRAD(1.0)
        IF(NJSET.LT.0) RETURN
C
C
C          Special initial state evolution for W-W fusion.
      ELSE
        CALL HEVOLV
        IF(NJSET.LT.0) RETURN
        DO 141 J=1,NJSET
141     JMATCH(J)=0
        DO 142 JET=1,2
          J=NJSET+1-2*JET
          PJSET(5,J)=-PJSET(5,JET)
142     JDCAY(J)=-2
        CALL QCDINI(NJSET-3,NJSET-1)
        IF(NJSET.LT.0) RETURN
      ENDIF
C
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      DO 200 I=3,NJSAVE,2
        JMATCH(I)=I+1
        JMATCH(I+1)=I
200   CONTINUE
      IF(NPAIR.EQ.0) THEN
        CALL DBLVEC(PJSET(1,N0JETS),DSUM)
        CALL DBLVEC(PJSET(1,N0JETS+1),DPASS)
        DO 231 K=1,4
231     DSUM(K)=DSUM(K)+DPASS(K)
        DSUM(5)=DSQRT(DSUM(4)**2-DSUM(1)**2-DSUM(2)**2-DSUM(3)**2)
        DO 232 K=1,5
          FRAME(K,1)=DSUM(K)
          FRAME(K,2)=FRAME(K,1)
232     CONTINUE
      ELSE
        DO 233 I=1,2
          DO 234 K=1,5
            FRAME(K,I)=PJSET(K,N0JETS+I-1)
234       CONTINUE
          IFRAME(I)=I
233     CONTINUE
      ENDIF
C
C          Set up and generate final state QCD parton shower.
C          Boost PJSET with -FRAME.
C
      DO 240 J=1,NJSAVE
        JET=JORIG(J)/JPACK
        IF(JET.EQ.0) JET=3
        IF(JET.GT.10) GO TO 240
C          Do this boost in double precision for 32-bit machines
        CALL DBOOST(-1,FRAME(1,JET),PJSET(1,J))
240   CONTINUE
C
C          Set maximum off-shell masses and JDCAY flags.
C
      IF(NPAIR.EQ.0) THEN
        NJFINL=N0JETS
        DO 340 J=N0JETS,NJSAVE
          IF(IABS(JTYPE(J)).LT.10) THEN
            PJSET(5,J)=EVOLMS(J,1.0)
            JDCAY(J)=-1
          ENDIF
340     CONTINUE
      ELSE
        NJFINL=N0PAIR
        DO 341 J=1,NPAIR
          IF(IABS(JTYPE(N0PAIR+J-1)).LT.10) THEN
            PJSET(5,N0PAIR+J-1)=PJETS(5,JPAIR(J))
            JDCAY(N0PAIR+J-1)=-1
          ENDIF
341     CONTINUE
      ENDIF
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END
