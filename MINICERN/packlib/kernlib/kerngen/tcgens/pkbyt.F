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
#include "sungs/pkbyt.F"
#elif !defined(CERNLIB_NOSHIFT)
      SUBROUTINE PKBYT (MIV,MBV,JTHP,NINTP,NBITS)
C
C CERN PROGLIB# M422    PKBYT           .VERSION KERNFOR  4.23  891215
C ORIG. 13/03/89  JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine
C
      DIMENSION    MIV(99), MBV(99), JTHP(9), NINTP(9), NBITS(2)

#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

      JTH  = JTHP(1)
      NINT = NINTP(1)
      IF (NINT.LE.0)         RETURN

      NZB  = NBITS(1)
      IF (NZB.GT.0)          GO TO 11
      NZB  = 1
      NPWD = NBITPW
      MSKA = 1
      GO TO 12

   11 NPWD = NBITS(2)
#if !defined(CERNLIB_QISASTD)
      MSKA = ISHFTR (IALL11, NBITPW-NZB)
#endif
#if defined(CERNLIB_QISASTD)
      MSKA = ISHFT  (IALL11,-NBITPW+NZB)
#endif

   12 JBV  = 1
      JIV  = 0
      IF (JTH.LT.2)          GO TO 21
      JBV  = (JTH-1)/NPWD + 1
      JPOS = JTH - (JBV-1)*NPWD - 1
      IF (JPOS.EQ.0)         GO TO 21
      NL   = JPOS*NZB
#if !defined(CERNLIB_QISASTD)
      MSKU = ISHFTL (MSKA,NL)
#endif
#if defined(CERNLIB_QISASTD)
      MSKU = ISHFT  (MSKA,NL)
#endif
      JIVE = NPWD - JPOS
      GO TO 22

C--                PACK EACH WORD

   21 NL   = 0
      MSKU = MSKA
      JIVE = JIV + NPWD
   22 JIVE = MIN (NINT,JIVE)
      IZW  = MBV(JBV)

   24 JIV  = JIV + 1
      IZW  = IOR (IAND(NOT(MSKU),IZW),
#if !defined(CERNLIB_QISASTD)
     +           IAND(MSKU,ISHFTL(MIV(JIV),NL)))
#endif
#if defined(CERNLIB_QISASTD)
     +           IAND(MSKU,ISHFT (MIV(JIV),NL)))
#endif
      IF (JIV.EQ.JIVE)       GO TO 27
      NL   = NL + NZB
#if !defined(CERNLIB_QISASTD)
      MSKU = ISHFTL (MSKU,NZB)
#endif
#if defined(CERNLIB_QISASTD)
      MSKU = ISHFT  (MSKU,NZB)
#endif
      GO TO 24

   27 MBV(JBV) = IZW
      IF (JIV.EQ.NINT)       RETURN
      JBV  = JBV + 1
      GO TO 21
      END
#endif
