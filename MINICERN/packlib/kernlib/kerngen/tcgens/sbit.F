*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:13  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/sbit.F"
#else
      SUBROUTINE SBIT (IT,IZW,IZP)
C
C CERN PROGLIB# M421    SBIT            .VERSION KERNFOR  4.23  891215
C MOD. true default 24/2/89, JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine

#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

#if !defined(CERNLIB_QISASTD)
      IZW = IOR ( IAND (IZW, NOT(ISHFTL(1,IZP-1)))
     +,           ISHFTL(IAND(IT,1), IZP-1) )
#endif
#if defined(CERNLIB_QISASTD)
      IZW = IOR ( IAND (IZW, NOT(ISHFT (1,IZP-1)))
     +,           ISHFT (IAND(IT,1), IZP-1) )
#endif
      RETURN
      END
#endif
