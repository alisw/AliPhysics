*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/02 23:17:49  thakulin
* Add support for EPC Fortran:  remove char-int equivalences and use
* F90 transfer facility instead.
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UHTOC (MS,NPW,MT,NCH)
C
C CERN PROGLIB# M409    UHTOC           .VERSION KERNFOR  4.21  890323
C ORIG. 10/02/89  JZ
C
      DIMENSION    MS(99)
      CHARACTER    MT*99

#include "kerngen/wordsize.inc"
      CHARACTER    CHWORD*(NCHAPW)
      INTEGER      IWORD
#if !defined(CERNLIB_F90) && !defined(CERNLIB_QFEPC)
      EQUIVALENCE (IWORD,CHWORD)
#endif

      IF   (NCH)             91, 29, 11
   11 IF (NPW.LE.0)          GO TO 91
      IF (NPW.EQ.1)          GO TO 21
      IF (NPW.LT.NCHAPW)     GO TO 31

C--------          NPW = maximum

      JT     = 0
      NWS    = NCH / NCHAPW
      NTRAIL = NCH - NWS*NCHAPW

      DO 14  JS=1,NWS
      IWORD  = MS(JS)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(JT+1:JT+NCHAPW) = transfer(IWORD,CHWORD)
#else
      MT(JT+1:JT+NCHAPW) = CHWORD
#endif
   14 JT  = JT + NCHAPW

      IF (NTRAIL.EQ.0)       RETURN

      IWORD = MS(NWS+1)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      CHWORD = transfer(IWORD,CHWORD)
#endif
      MT(JT+1:JT+NTRAIL) = CHWORD(1:NTRAIL)
      RETURN

C--------          NPW = 1
C--                equivalent to 'CALL UH1TOC(MS,MT,NCH)'

   21 DO 24  JS=1,NCH
      IWORD  = MS(JS)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      CHWORD = transfer(IWORD,CHWORD) 
#endif
      MT(JS:JS) = CHWORD(1:1)
   24 CONTINUE
   29 RETURN

C--------          NPW = 2 ...

   31 JT     = 0
      NWS    = NCH / NPW
      NTRAIL = NCH - NWS*NPW

      DO 34  JS=1,NWS
      IWORD  = MS(JS)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      CHWORD = transfer(IWORD,CHWORD)
#endif
      MT(JT+1:JT+NPW) = CHWORD(1:NPW)
   34 JT  = JT + NPW

      IF (NTRAIL.EQ.0)       RETURN

      IWORD = MS(NWS+1)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      CHWORD = transfer(IWORD,CHWORD)
#endif
      MT(JT+1:JT+NTRAIL) = CHWORD(1:NTRAIL)
      RETURN

   91 CALL ABEND
      END
