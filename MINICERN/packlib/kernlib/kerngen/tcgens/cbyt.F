*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:10  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/cbyt.F"
#elif defined(CERNLIB_QMLNX)
#include "lnxgs/cbyt.F"
#else
      SUBROUTINE CBYT (M,LM,IX,LX,NB)
C
C CERN PROGLIB# M421    CBYT            .VERSION KERNFOR  4.23  891215
C ORIG. 14/08/70 CL
C
#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_jbit.inc"
* Ignoring t=pass
#endif
C
#if !defined(CERNLIB_QISASTD)
      CALL SBYT (JBYT(M,LM,NB),IX,LX,NB)
#endif
#if defined(CERNLIB_QISASTD)
      IX = IOR ( IAND (IX,
     +                  NOT(ISHFT (ISHFT(IALL11,-NBITPW+NB), LX-1)))
     +,   ISHFT (ISHFT(ISHFT(M,-LM+1),NBITPW-NB), -NBITPW-1+LX+NB))
#endif
      RETURN
      END
#endif
