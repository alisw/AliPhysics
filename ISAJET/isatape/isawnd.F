#include "isajet/pilot.h"
      SUBROUTINE ISAWND
C
C          WRITE END RECORD, TYPE 300
C          CONTAINS CROSS SECTIONS AND LUMINOSITY
C
#include "isajet/itapes.inc"
#include "isajet/final.inc"
#include "isajet/totals.inc"
#include "isajet/zevel.inc"
#include "isajet/jetlim.inc"
#include "isajet/dylim.inc"
#include "isajet/keys.inc"
C
C       keep entry point WREND for backward compatibility
      ENTRY WREND
      ITA=IABS(ITEVT)
      IZEVEL(1)=300
      IZEVEL(2)=1
      IZEVEL(3)=NKINF
      ZEVEL(4)=SIGF
      ZEVEL(5)=ALUM
      ZEVEL(6)=ACCEPT
      IZEVEL(7)=NRECS
      IL=7
      CALL BUFOUT(IL)
      WRITE(ITLIS,1010) NRECS,ITA
1010  FORMAT(////' THIS RUN WROTE',I10,
     1' PHYSICAL RECORDS ON TAPE',I3)
      RETURN
      END
