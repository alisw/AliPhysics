*
* $Id$
*
* $Log$
* Revision 1.4  1998/09/25 09:32:02  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.3  1997/09/02 14:27:00  mclareni
* WINNT correction
*
* Revision 1.2  1997/02/04 17:36:25  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:31:40  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QFMSOFT)
#include "wntgs/ublow.F"
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#include "dosgs/ublow.F"
#elif defined(CERNLIB_QMMPW)
#include "mpwgs/ublow.F"
#elif defined(CERNLIB_QMVAOS)||defined(CERNLIB_QMVMI)
#include "allgs/ublow.F"
#elif (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
#include "lnxgs/ublow.F"
#elif (defined(CERNLIB_QMLNX) && defined(CERNLIB_PPC))
#include "lnxppcgs/ublow.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/ublow.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/ublow.F"
#elif defined(CERNLIB_B32)||defined(CERNLIB_B64)
      SUBROUTINE UBLOW (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBLOW           .VERSION KERNFOR  4.30  910819
C ORIG. 05/12/89, FCA+JZ
C

      DIMENSION    MS(99), MT(99), NCHP(9)
#include "kerngen/iallbl.inc"
#include "kerngen/ublowx1.inc"

      NCH = NCHP(1)
      IF   (NCH)             91, 29, 11
#if defined(CERNLIB_B64)
   11 NWS    = ISHFT (NCH,-3)
      NTRAIL = IAND (NCH,7)
#endif
#if defined(CERNLIB_B32)
   11 NWS    = ISHFT (NCH,-2)
      NTRAIL = IAND (NCH,3)
#endif
      JT     = 0
      IF (NWS.EQ.0)          GO TO 26

C--                Unpack the initial complete words

      DO 24 JS=1,NWS
      MWD      = MS(JS)
      MT(JT+1) = IOR (IBLAN1,IAND(MASK1,MWD))
      MT(JT+2) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD, 8)))
      MT(JT+3) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD,16)))
#if defined(CERNLIB_B64)
      MT(JT+4) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD,24)))
      MT(JT+5) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD,32)))
      MT(JT+6) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD,40)))
      MT(JT+7) = IOR (IBLAN1,IAND(MASK1,ISHFT(MWD,48)))
      MT(JT+8) = IOR (IBLAN1,           ISHFT(MWD,56) )
   24 JT = JT + 8
#endif
#if defined(CERNLIB_B32)
      MT(JT+4) = IOR (IBLAN1,           ISHFT(MWD,24) )
   24 JT = JT + 4
#endif

      IF (NTRAIL.EQ.0)       RETURN

C--                Unpack the trailing word

   26 MWD = MS(NWS+1)

      DO 28 JS=1,NTRAIL
      MT(JT+1) = IOR (IBLAN1,IAND(MASK1,MWD))
      MWD = ISHFT (MWD,8)
   28 JT = JT + 1
   29 RETURN

   91 CALL ABEND
      END
#endif
