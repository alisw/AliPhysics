#include "isajet/pilot.h"
      SUBROUTINE IDGEN
C
C          Call system date and time routines (non-standard) to set up
C          run identification:
C          IDVER=100*VERSN     (integer ISAJET version number)
C          IDG(1)=YYMMDD       (integer year-month-day)
C          IDG(2)=HHMMSS       (integer hour-minute-second)
C
#include "isajet/itapes.inc"
#include "isajet/idrun.inc"
#if defined(CERNLIB_CDC)||defined(CERNLIB_ETA)
      CHARACTER*10 CHAR,DATE,TIME
#endif
#if defined(CERNLIB_SUN)||defined(CERNLIB_SGI)
      DIMENSION ISUN(3)
#endif
C          Default run id is zero.
      IYMD=0.
      IHMS=0.
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_NOCERN))
C          Call CDC date and time and convert to integer.
      CHAR=DATE()
      READ(CHAR,'(1X,I2,1X,I2,1X,I2,1X)') IA,IB,IC
      IYMD=10000*IC+100*IA+IB
      CHAR=TIME()
      READ(CHAR,'(1X,I2,1X,I2,1X,I2,1X)') IA,IB,IC
      IHMS=10000*IA+100*IB+IC
#endif
#if (defined(CERNLIB_ETA))&&(defined(CERNLIB_NOCERN))
C          Call ETA date and time and convert to integer.
      CHAR=DATE()
      READ(CHAR,'(I2,1X,I2,1X,I2)') IA,IB,IC
      IYMD=10000*IC+100*IA+IB
      CHAR=TIME()
      READ(CHAR,'(I2,1X,I2,1X,I2)') IA,IB,IC
      IHMS=10000*IA+100*IB+IC
#endif
#if (defined(CERNLIB_SGI))&&(defined(CERNLIB_NOCERN))
C          Call Silicon Graphics date and time
      CALL IDATE(ISUN(1),ISUN(2),ISUN(3))
      IYMD=10000*ISUN(3)+100*ISUN(2)+ISUN(1)
      CALL ITIME(ISUN)
      IHMS=10000*ISUN(1)+100*ISUN(2)+ISUN(3)
#endif
#if (defined(CERNLIB_SUN))&&(defined(CERNLIB_NOCERN))
C          Call SUN date and time
      CALL IDATE(ISUN)
      IYMD=10000*(MOD(ISUN(3),100))+100*ISUN(2)+ISUN(1)
      CALL ITIME(ISUN)
      IHMS=10000*ISUN(1)+100*ISUN(2)+ISUN(3)
#endif
#if defined(CERNLIB_IBM)||defined(CERNLIB_VAX)||defined(CERNLIB_CERN)
C          Call DATIME for date and time. (In Cern library)
      CALL DATIME(IYMD,IHMS)
#endif
      IDG(1)=IYMD
      IDG(2)=IHMS
      RETURN
      END
