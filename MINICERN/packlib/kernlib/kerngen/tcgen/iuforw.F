*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/iuforw.F"
#else
      FUNCTION IUFORW (IC,JL,JR)
C
C CERN PROGLIB# M410    IUFORW          .VERSION KERNFOR  4.23  891215
C ORIG. 27/05/67
C
      COMMON /SLATE/ ND,JEND, DUMMY(38)
      DIMENSION    IC(80)
#if !defined(CERNLIB_QASCII)||!defined(CERNLIB_QISASTD)||defined(CERNLIB_B36)||defined(CERNLIB_B48)||defined(CERNLIB_B60)
      DIMENSION    NUM(11)
      DATA  NUM   /4H    ,4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,4H5   ,
     +                           4H6   ,4H7   ,4H8   ,4H9    /
#endif

      ND   = 0
      JEND = JL
      IVAL = 0
C
   12 IF (JEND.GT.JR)        GO TO 20
#if !defined(CERNLIB_QASCII)||!defined(CERNLIB_QISASTD)||defined(CERNLIB_B36)||defined(CERNLIB_B48)||defined(CERNLIB_B60)
      K = IUCOMP (IC(JEND),NUM,11) - 2
      IF (K.EQ.-1)           GO TO 19
#endif
#if (defined(CERNLIB_QASCII))&&(defined(CERNLIB_QISASTD))&&(defined(CERNLIB_QORTHOLL))&&(defined(CERNLIB_B64))
      K = ISHFT (IC(JEND),-56)
#endif
#if (defined(CERNLIB_QASCII))&&(defined(CERNLIB_QISASTD))&&(defined(CERNLIB_QORTHOLL))&&(defined(CERNLIB_B32))
      K = ISHFT (IC(JEND),-24)
#endif
#if (defined(CERNLIB_QASCII))&&(defined(CERNLIB_QISASTD))&&(!defined(CERNLIB_QORTHOLL))&&(defined(CERNLIB_B32)||defined(CERNLIB_B64))
      K = IAND (IC(JEND),255)
#endif
#if (defined(CERNLIB_QASCII))&&(defined(CERNLIB_QISASTD))&&(defined(CERNLIB_B32)||defined(CERNLIB_B64))
      IF (K.EQ.32)           GO TO 19
      IF (K.GE.58)           GO TO 20
      K = K - 48
#endif
      IF (K.LT.0)            GO TO 20
      IVAL = 10*IVAL + K
      ND   = ND + 1
   19 JEND = JEND + 1
      GO TO 12
C
   20 IUFORW = IVAL
      RETURN
      END
#endif
