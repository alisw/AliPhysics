*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QIEEE)
      FUNCTION IUWEED (IA,N)
C
C CERN PROGLIB# N103    IUWEED          .VERSION KERNFOR  4.31  911111
C ORIG. 16/10/91 JZ
C-    Find position of first number which is infinity or not-a-number
C-    in the N word vector IA, return zero if none
C-    ---  this version only for machines with IEEE representation  ---
      DIMENSION IA(*)
#if defined(CERNLIB_B32)
      PARAMETER (NMANT=23, IVBAD=255)
#endif
#if defined(CERNLIB_B64)
      PARAMETER (NMANT=52, IVBAD=2047)
#endif
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
#endif

      DO  24  J=1,N
#if !defined(CERNLIB_QISASTD)
      IVAL = IAND (ISHFTR(IA(J),NMANT), IVBAD)
#endif
#if defined(CERNLIB_QISASTD)
      IVAL = IAND (ISHFT(IA(J),-NMANT), IVBAD)
#endif
      IF (IVAL.EQ.IVBAD)           GO TO 29
   24 CONTINUE
      J = 0
   29 IUWEED = J
      RETURN
      END
#endif
