#include "isajet/pilot.h"
      SUBROUTINE EVOL01
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 1 (TWOJET)
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
#include "isajet/frame.inc"
      REAL    EVOLMS
      INTEGER I,K,J,NJSAVE,NJFINL
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
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
      CALL IFRAMS(N0JETS,NJSAVE,1,.FALSE.)
C
C          Set maximum off-shell masses and JDCAY flags.
C
      NJFINL=N0JETS
      DO 310 J=N0JETS,NJSAVE
        IF(IABS(JTYPE(J)).LT.10) THEN
          PJSET(5,J)=EVOLMS(J,1.0)
          JDCAY(J)=-1
        ENDIF
310   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END
