*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/jbytet.F"
#else
      FUNCTION JBYTET (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    JBYTET          .VERSION KERNFOR  4.23  891215
C ORIG. 13/03/89  JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine
#include "kerngen/wordsize.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

#if !defined(CERNLIB_QISASTD)
      JBYTET = IAND (MZ, ISHFTR(ISHFTL(IZW,NBITPW+1-IZP-NZB),
     +                           NBITPW-NZB) )
#endif
#if (defined(CERNLIB_QISASTD))&&(!defined(CERNLIB_QMILSTD))
      JBYTET = IAND (MZ, ISHFT (ISHFT (IZW,NBITPW+1-IZP-NZB),
     +                          -NBITPW+NZB) )
#endif
#if defined(CERNLIB_QMILSTD)
      JBYTET = IAND (MZ, IBITS(IZW,IZP-1,NZB))
#endif
      RETURN
      END
#endif
