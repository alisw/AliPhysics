*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:14  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/sbytor.F"
#else
      SUBROUTINE SBYTOR (IT,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    SBYTOR          .VERSION KERNFOR  4.23  891215
C ORIG. 13/03/89  JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine
C
#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

#if !defined(CERNLIB_QISASTD)
      IZW = IOR (IZW,
     +      ISHFTL (IAND(IT,ISHFTR(IALL11, NBITPW-NZB)), IZP-1))
#endif
#if defined(CERNLIB_QISASTD)
      IZW = IOR (IZW,
     +      ISHFT  (IAND(IT,ISHFT (IALL11,-NBITPW+NZB)), IZP-1))
#endif
      RETURN
      END
#endif
