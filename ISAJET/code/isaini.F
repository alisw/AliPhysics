#include "isajet/pilot.h"
      SUBROUTINE ISAINI(JTDKY,JTEVT,JTCOM,JTLIS)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-         INITIALIZE PROCESSES
C-
C-   Inputs  : 
C          JTDKY = +/- UNIT NUMBER FOR DECAY TABLE FILE.
C                      IF IT IS NEGATIVE, DECAY TABLE IS NOT PRINTED.
C          JTEVT = +/- UNIT NUMBER FOR OUTPUT EVENT FILE.
C                      IF IT IS NEGATIVE, ONLY STABLE PARTICLES ARE
C                      WRITTEN ON IT.
C          JTCOM =     UNIT NUMBER FOR COMMAND FILE.
C          JTLIS =     UNIT NUMBER FOR LISTING.
C-
C-   Created   3-FEB-1988   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
C
#include "isajet/idrun.inc"
#include "isajet/itapes.inc"
C
C          ENTRY.
      ITDKY=IABS(JTDKY)
      ITEVT=JTEVT
      ITCOM=IABS(JTCOM)
      ITLIS=IABS(JTLIS)
C
      IEVT=0
      CALL SETCON
      CALL RESET
      IF(JTDKY.GT.0) THEN
        CALL SETDKY(.TRUE.)
      ELSE
        CALL SETDKY(.FALSE.)
      ENDIF
C
  999 RETURN
      END
