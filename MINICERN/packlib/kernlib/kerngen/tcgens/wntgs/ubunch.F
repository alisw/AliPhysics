*
* $Id$
*
* $Log$
* Revision 1.2  1997/09/02 14:27:07  mclareni
* WINNT correction
*
* Revision 1.1  1997/02/04 17:36:30  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1  1996/02/15 17:50:24  mclareni
* Kernlib
*
*
#if !defined(CERNLIB_QF2C)
      SUBROUTINE UBUNCH (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBUNCH          .VERSION KERNDOS  1.00  920624
C ORIG. 03/02/89 K.M.STORR
C

      DIMENSION    MS(99), MT(99), NCHP(9)
#ifndef CERNLIB_QFMSOFT
      PARAMETER    (IBLANK = X'20202020')
      PARAMETER    (MASK1  = X'000000FF')
#else
      PARAMETER    (IBLANK = 16#20202020)
      PARAMETER    (MASK1  = 16#000000FF)
#endif

#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"

      NCH = NCHP(1)
      IF   (NCH)             91,39,11
   11 NWT    = ishft (NCH,-2)
      NTRAIL = IAND (NCH,3)
      JS     = 0
      IF (NWT.EQ.0)          GO TO 31

C--                Pack the initial complete words

      DO 24  JT=1,NWT
#ifndef CERNLIB_QF2C
      MT(JT) = IOR (IOR (IOR (
     +                  IAND(MS(JS+1),MASK1),
     +           ishft (IAND(MS(JS+2),MASK1), 8)),
     +           ishft (IAND(MS(JS+3),MASK1),16)),
     +           ishft      (MS(JS+4),       24) )
#else
      MT(JT) = OR (OR (OR (
     +                  AND(MS(JS+1),MASK1),
     +           LshIft (AND(MS(JS+2),MASK1), 8)),
     +           LshIft (AND(MS(JS+3),MASK1),16)),
     +           LshIft     (MS(JS+4),       24) )
#endif
   24 JS = JS + 4

      IF (NTRAIL.EQ.0)       RETURN

C--                Pack the trailing word

   31 MWD = IBLANK
      JS  = NCH

      DO 34 JT=1,NTRAIL
#ifndef CERNLIB_QF_F2C
      MWD = IOR (ishft(MWD,8), IAND(MS(JS),MASK1))
#else
      MWD = OR (LshIft(MWD,8), AND(MS(JS),MASK1))
#endif
   34 JS  = JS - 1
      MT(NWT+1) = MWD
   39 RETURN

   91 CALL ABEND
      END
#endif
