#include "isajet/pilot.h"
      SUBROUTINE ISTRAD(FUDGE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      Set parameters and call QCDINI to generate initial
C-      state radiation
C-   Inputs  : 
C-     FUDGE= fudge factor
C-
C-   Created  16-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL    FUDGE
#include "isajet/limevl.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/jetpar.inc"
      REAL    OFF
      INTEGER I
C----------------------------------------------------------------------
C
      IF ( USELIM.AND.CONCUT.LT.1.0 ) THEN
        OFF=ETTHRS
      ELSEIF( .NOT.USELIM) THEN
        OFF=SQRT(QSQ)*FUDGE
      ELSE
        OFF=SQRT(QSQ)
      ENDIF
      DO 150 I=1,2
        PJSET(5,I)=-OFF
150   JDCAY(I)=-2
      JMATCH(1)=0
      JMATCH(2)=0
C
      CALL QCDINI(1,2)
  999 RETURN
      END
