*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSETHI (INTP, CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CSETHI          .VERSION KERNFOR  4.31  911111
C ORIG. 17/10/89, JZ
C
C-    Set hexadecimal integer into CHV(JL:JR) right-justified

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
      K    = IAND (IVAL,15)
#if !defined(CERNLIB_QISASTD)
      IVAL = ISHFTR (IVAL,4)
#endif
#if defined(CERNLIB_QISASTD)
      IVAL = ISHFT (IVAL,-4)
#endif
      IF (K.LT.10)  THEN
#if defined(CERNLIB_QASCII)
          CHV(JJ:JJ) = CHAR(K+48)
        ELSE
          CHV(JJ:JJ) = CHAR(K+55)
#endif
#if defined(CERNLIB_QEBCDIC)
          CHV(JJ:JJ) = CHAR(K+240)
        ELSE
          CHV(JJ:JJ) = CHAR(K+183)
#endif
        ENDIF

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
