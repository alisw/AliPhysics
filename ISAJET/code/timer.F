#include "isajet/pilot.h"
      SUBROUTINE TIMER(IT)
C
C          CALL SYSTEM CPU CLOCK -- MACHINE DEPENDENT.
C          IT=1 FOR RUN START TIME.
C          IT=2 FOR RUN STOP TIME.
C
#include "isajet/itapes.inc"
#include "isajet/times.inc"
      DIMENSION TIMES(2)
      EQUIVALENCE (TIMES(1),TIME1)
      DIMENSION TTT(2)
#if defined(CERNLIB_VAX)
      INTEGER CPUTIM(2),ITMLST(4),NHSEC
      EXTERNAL JPI$_CPUTIM
#endif
C
C          DEFAULT IS TO RETURN ZERO.
      TNOW=0.
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_NOCERN))
C          SECOND GIVES CPU TIME ON CDC.
      CALL SECOND(TNOW)
#endif
#if (defined(CERNLIB_ETA))&&(defined(CERNLIB_NOCERN))
C          SECOND GIVES CPU TIME ON ETA.
      TNOW=SECOND()
#endif
#if (defined(CERNLIB_IBMRT))&&(defined(CERNLIB_NOCERN))
C          MCLOCK GIVES CPU TIME ON IBM RS/6000.
      TNOW=FLOAT(MCLOCK())/60.
#endif
#if (defined(CERNLIB_SGI))&&(defined(CERNLIB_NOCERN))
C          ETIME GIVES CPU TIME ON SILICON GRAPHICS.
      TNOW=ETIME(TTT)
#endif
#if (defined(CERNLIB_SUN))&&(defined(CERNLIB_NOCERN))
C          ETIME GIVES CPU TIME ON SUN.
      TNOW=ETIME(TTT)
#endif
#if (defined(CERNLIB_VAX))&&(defined(CERNLIB_NOCERN))
C          VAX HAS NO FORTRAN FUNCTION FOR CPU TIME.
C          FOLLOWING PROVIDED BY T. KILLIAN
      ITMLST(1)=ISHFT(%LOC(JPI$_CPUTIM),16)+4
      ITMLST(2)=%LOC(NHSEC)
      ITMLST(3)=0
      ITMLST(4)=0
      CALL SYS$GETJPI(,,,ITMLST,,,)
      TNOW=.01*NHSEC
#endif
#if defined(CERNLIB_CERN)
      CALL TIMEST(1.E7)
      CALL TIMEX(TNOW)
#endif
      TIMES(IT)=TNOW
      RETURN
      END
