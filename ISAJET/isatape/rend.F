#include "isajet/pilot.h"
      SUBROUTINE REND
C
C          INVERSE OF WREND
C          READ END RECORD (TYPE 300)
C
#include "isajet/itapes.inc"
#include "isajet/final.inc"
#include "isajet/zevel.inc"
      NKINF=IZEVEL(3)
      SIGF=ZEVEL(4)
      ALUM=ZEVEL(5)
      ACCEPT=ZEVEL(6)
      NRECS=IZEVEL(7)
      RETURN
      END
