*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSETOI (INTP, CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CSETOI          .VERSION KERNFOR  4.33  920511
C ORIG. 03/06/92, JZ
C
C-    Set octal integer into CHV(JL:JR) right-justified

      DIMENSION    INTP(9), JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)

#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
#endif

      JL  = JLP(1)
      JJ  = JRP(1)

      IVAL   = INTP(1)
      NDG    = 0
      NGSLAT = 0

   12 IF (JJ.LT.JL)          GO TO 97
      K    = IAND (IVAL,7)
#if !defined(CERNLIB_QISASTD)
      IVAL = ISHFTR (IVAL,3)
#endif
#if defined(CERNLIB_QISASTD)
      IVAL = ISHFT (IVAL,-3)
#endif
#if defined(CERNLIB_QASCII)
      CHV(JJ:JJ) = CHAR(K+48)
#endif
#if defined(CERNLIB_QEBCDIC)
      CHV(JJ:JJ) = CHAR(K+240)
#endif

      JJ  = JJ  - 1
      NDG = NDG + 1
      IF (IVAL.NE.0)         GO TO 12
      GO TO 98

   97 NGSLAT = JL
   98 NFSLAT = JJ
      NESLAT = JJ
      NDSLAT = NDG
      RETURN
      END
