#include "isajet/pilot.h"
      SUBROUTINE ISABEG(IFL)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       Initialize a process before event generation
C-
C-   Created   5-FEB-1988   Serban D. Protopopescu
C-
C    Ver 7.14: Do logic after setting physics parameters
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/nodcay.inc"
#include "isajet/idrun.inc"
#include "isajet/keys.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/isloop.inc"
#include "isajet/xmssm.inc"
#include "isajet/isapw.inc"
C
      INTEGER IFL,I
      LOGICAL FIRST
      SAVE FIRST
      CHARACTER*30 ISAPW2
      SAVE ISAPW2
      DATA FIRST/.TRUE./
C          ISAPW2 is used to check whether ALDATA is loaded
      DATA ISAPW2/'ALDATA REQUIRED BY FORTRAN G,H'/
C
C          Initialize
C
      IF(ISAPW1.NE.ISAPW2) THEN
        PRINT*, ' ISABEG ERROR: BLOCK DATA ALDATA HAS NOT BEEN LOADED.'
        PRINT*, ' ISAJET CANNOT RUN WITHOUT IT.'
        PRINT*, ' PLEASE READ THE FINE MANUAL FOR ISAJET.'
        STOP99
      ENDIF
C
      IF (FIRST) THEN
        FIRST=.FALSE.
      ELSE
        CALL SETNXT
      ENDIF
      IEVT=0
      IEVGEN=0
      NEVENT=0
      IEVOL=1
      IFRG=1
C
C          Read in user data and decay table
C

      CALL READIN(IFL)
      IF(IFL.NE.0) GOTO 999
      CALL IDGEN
      IF(GOMSSM) THEN
        CALL DOMSSM
      ENDIF
      IF (KEYS(10).AND..NOT.GOMSSM) THEN
        CALL SETH
      END IF
      CALL SETDKY(.FALSE.)
C
C          Generate NSIGMA unevolved events for SIGF calculation
C
C          TWOJET events
      IF(KEYS(1)) THEN
        CALL MBSET
        CALL SETW
        CALL LOGIC
        CALL PRTLIM
        CALL PTFUN
        DO 105 I=1,NSIGMA
105     CALL TWOJET
        CALL TIMER(1)
C
C          E+E- events
      ELSE IF(KEYS(2)) THEN
        CALL SETW
        CALL LOGIC
        CALL PRTLIM
        CALL EEBEG
        CALL EEMAX
        DO 205 I=1,NSIGMA
205     CALL ELCTRN
        CALL TIMER(1)
C
C          DRELLYAN events
      ELSE IF(KEYS(3)) THEN
        CALL SETW
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL QFUNC
        DO 305 I=1,NSIGMA
305     CALL DRLLYN
        CALL TIMER(1)
C
C          MINBIAS events
      ELSE IF(KEYS(4)) THEN
        PBEAM(1)=HALFE
        PBEAM(2)=HALFE
        CALL PRTLIM
        CALL MBSET
        CALL TIMER(1)
C
C          SUPERSYM events
      ELSE IF(KEYS(5)) THEN
        CALL SETW
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL PTFUN
        DO 505 I=1,NSIGMA
505     CALL TWOJET
        CALL TIMER(1)
C
C          WPAIR events
      ELSE IF(KEYS(6)) THEN
        CALL SETW
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL PTFUN
        DO 605 I=1,NSIGMA
        CALL TWOJET
605     CALL WPAIR
        CALL TIMER(1)
C
C          HIGGS events
      ELSE IF(KEYS(7)) THEN
        CALL SETW
        IF(GOMSSM) THEN
          CALL SETHSS
        ELSE
          CALL SETH
        ENDIF
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL QFUNC
        DO 705 I=1,NSIGMA
705     CALL DRLLYN
        CALL TIMER(1)
C
C          PHOTON events
      ELSEIF(KEYS(8)) THEN
        CALL MBSET
        CALL SETW
        CALL LOGIC
        CALL PRTLIM
        CALL PTFUN
        DO 805 I=1,NSIGMA
805     CALL TWOJET
        CALL TIMER(1)
C
C          TCOLOR events
      ELSE IF(KEYS(9)) THEN
        CALL SETW
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL QFUNC
        DO 905 I=1,NSIGMA
905     CALL DRLLYN
        CALL TIMER(1)
C
C          WHIGGS events
      ELSE IF(KEYS(10)) THEN
        CALL SETW
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL PTFUN
        DO 906 I=1,NSIGMA
        CALL TWOJET
906     CALL WHIGGS
        CALL TIMER(1)
C
C          EXTRADIM events
      ELSE IF(KEYS(11)) THEN
        CALL SETW
        CALL SETKKG
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL QFUNC
        DO 1105 I=1,NSIGMA
          CALL DRLLYN
1105    CONTINUE
        CALL TIMER(1)
C
C          ZJJ events
C          ZJJ0 initializes cross sections, so no event loop
      ELSEIF(KEYS(12)) THEN
        CALL SETW
        CALL MGINIT
        CALL MBSET
        CALL LOGIC
        CALL PRTLIM
        CALL ZJJ0
        CALL TIMER(1)
      ELSE
        STOP 99
      ENDIF
999   RETURN
      END
