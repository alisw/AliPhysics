*
* $Id$
*
* $Log$
* Revision 1.4  1998/09/25 09:32:06  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.3  1997/09/02 14:27:00  mclareni
* WINNT correction
*
* Revision 1.2  1997/02/04 17:36:25  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:31:41  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QFMSOFT) && defined(CERNLIB_WINNT)
#include "wntgs/ubunch.F"
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#include "dosgs/ubunch.F"
#elif defined(CERNLIB_QMVAOS)||defined(CERNLIB_QMVMI)
#include "allgs/ubunch.F"
#elif (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
#include "lnxgs/ubunch.F"
#elif (defined(CERNLIB_QMLNX) && defined(CERNLIB_PPC))
#include "lnxppcgs/ubunch.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/ubunch.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/ubunch.F"
#elif defined(CERNLIB_B32)||defined(CERNLIB_B64)
      SUBROUTINE UBUNCH (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBUNCH          .VERSION KERNFOR  4.30  910819
C ORIG. 05/12/89, FCA+JZ
C

      DIMENSION    MS(99), MT(99), NCHP(9)
#include "kerngen/iallbl.inc"
#include "kerngen/ubnchx1.inc"

      NCH = NCHP(1)
      IF   (NCH)             91,39,11
#if defined(CERNLIB_B64)
   11 NWT    = ISHFT (NCH,-3)
      NTRAIL = IAND (NCH,7)
#endif
#if defined(CERNLIB_B32)
   11 NWT    = ISHFT (NCH,-2)
      NTRAIL = IAND (NCH,3)
#endif
      JS     = 0
      IF (NWT.EQ.0)          GO TO 31

C--                Pack the initial complete words

#if defined(CERNLIB_B64)
      DO 24  JT=1,NWT
      MT(JT) = IOR (IOR (IOR (IOR (IOR (IOR (IOR (
     +                  IAND(MS(JS+1),MASK1),
     +           ISHFT (IAND(MS(JS+2),MASK1), -8)),
     +           ISHFT (IAND(MS(JS+3),MASK1),-16)),
     +           ISHFT (IAND(MS(JS+4),MASK1),-24)),
     +           ISHFT (IAND(MS(JS+5),MASK1),-32)),
     +           ISHFT (IAND(MS(JS+6),MASK1),-40)),
     +           ISHFT (IAND(MS(JS+7),MASK1),-48)),
     +           ISHFT      (MS(JS+8),       -56) )
   24 JS = JS + 8
#endif
#if defined(CERNLIB_B32)
      DO 24  JT=1,NWT
      MT(JT) = IOR (IOR (IOR (
     +                  IAND(MS(JS+1),MASK1),
     +           ISHFT (IAND(MS(JS+2),MASK1), -8)),
     +           ISHFT (IAND(MS(JS+3),MASK1),-16)),
     +           ISHFT      (MS(JS+4),       -24) )
   24 JS = JS + 4
#endif

      IF (NTRAIL.EQ.0)       RETURN

C--                Pack the trailing word

   31 MWD = IALLBL
      JS  = NCH

      DO 34 JT=1,NTRAIL
      MWD = IOR (ISHFT(MWD,-8), IAND(MS(JS),MASK1))
   34 JS  = JS - 1
      MT(NWT+1) = MWD
   39 RETURN

   91 CALL ABEND
      END
#endif
