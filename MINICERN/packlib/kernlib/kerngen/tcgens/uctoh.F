*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/02 23:17:45  thakulin
* Add support for EPC Fortran:  remove char-int equivalences and use
* F90 transfer facility instead.
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UCTOH (MS,MT,NPW,NCH)
C
C CERN PROGLIB# M409    UCTOH           .VERSION KERNFOR  4.21  890323
C ORIG. 10/02/89  JZ
C
      CHARACTER    MS*99
      DIMENSION    MT(99)

#include "kerngen/wordsize.inc"
      CHARACTER    CHWORD*(NCHAPW)
      CHARACTER    BLANK *(NCHAPW)
      PARAMETER   (BLANK = ' ')
      INTEGER      IWORD
#if !defined(CERNLIB_F90) && !defined(CERNLIB_QFEPC)
      EQUIVALENCE (IWORD,CHWORD)
#endif

      IF   (NCH)             91, 29, 11
   11 IF (NPW.LE.0)          GO TO 91
      IF (NPW.EQ.1)          GO TO 21
      IF (NPW.LT.NCHAPW)     GO TO 31

C--------          NPW = maximum

      JS     = 0
      NWT    = NCH / NCHAPW
      NTRAIL = NCH - NWT*NCHAPW

      DO 14  JT=1,NWT
      CHWORD = MS(JS+1:JS+NCHAPW)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(JT) = transfer(CHWORD,IWORD)
#else
      MT(JT) = IWORD
#endif
   14 JS     = JS + NCHAPW

      IF (NTRAIL.EQ.0)       RETURN

      CHWORD    = MS(JS+1:JS+NTRAIL)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(NWT+1) = transfer(CHWORD,IWORD)
#else
      MT(NWT+1) = IWORD
#endif
      RETURN

C--------          NPW = 1
C--                equivalent to 'CALL UCTOH1(MS,MT,NCH)'

   21 CHWORD = BLANK

      DO 24 JS=1,NCH
      CHWORD(1:1) = MS(JS:JS)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(JS)      = transfer(CHWORD,IWORD)
#else
      MT(JS)      = IWORD
#endif
   24 CONTINUE
   29 RETURN

C--------          NPW = 2 ...

   31 CHWORD = BLANK

      JS     = 0
      NWT    = NCH / NPW
      NTRAIL = NCH - NWT*NPW

      DO 34  JT=1,NWT
      CHWORD(1:NPW) = MS(JS+1:JS+NPW)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(JT) = transfer(CHWORD,IWORD)
#else
      MT(JT) = IWORD
#endif
   34 JS     = JS + NPW

      IF (NTRAIL.EQ.0)       RETURN

      CHWORD    = MS(JS+1:JS+NTRAIL)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(NWT+1) = transfer(CHWORD,IWORD)
#else
      MT(NWT+1) = IWORD
#endif
      RETURN

   91 CALL ABEND
      END
