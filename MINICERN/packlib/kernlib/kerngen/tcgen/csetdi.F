*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:43  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CSETDI (INTP, CHV,JLP,JRP)
C
C CERN PROGLIB# M432    CSETDI          .VERSION KERNFOR  4.29  910718
C ORIG. 04/10/88, JZ
C
C-    Set decimal integer into CHV(JL:JR) right-justified

      DIMENSION    INTP(9), JLP(9), JRP(9)

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHV*(*)


      INT = INTP(1)
      JL  = JLP(1)
      JJ  = JRP(1)

      IVAL   = IABS (INT)
      NDG    = 0
      NGSLAT = 0

   12 IF (JJ.LT.JL)          GO TO 97
      K    = MOD (IVAL,10)
      IVAL = IVAL/10
#if defined(CERNLIB_QASCII)
      CHV(JJ:JJ) = CHAR(K+48)
#endif
#if defined(CERNLIB_QEBCDIC)
      CHV(JJ:JJ) = CHAR(K+240)
#endif

      JJ  = JJ  - 1
      NDG = NDG + 1
      IF (IVAL.NE.0)         GO TO 12
      IF (INT.GE.0)          GO TO 98
      IF (JJ.LT.JL)          GO TO 97
      CHV(JJ:JJ) = '-'
      NFSLAT = JJ - 1
      GO TO 99

   97 NGSLAT = JL
   98 NFSLAT = JJ
   99 NESLAT = JJ
      NDSLAT = NDG
      RETURN
      END
