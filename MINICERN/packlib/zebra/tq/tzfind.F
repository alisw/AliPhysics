*
* $Id$
*
* $Log$
* Revision 1.3  1997/03/14 17:22:21  mclareni
* WNT mods
*
* Revision 1.2.2.1  1997/01/21 11:34:12  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.2  1996/04/18 16:14:51  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZFIND (IXSTOR,LBK,IDHP,IDNP,IFLAG)

C-    Find title bank with IDH (and IDN if non-zero)

#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
C--------------    END CDE                             -----------------
      INTEGER      IDHP(9), IDNP(9), IDH
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HTZFI, 4HND   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HTZFIND /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'TZFIND  ')
#endif

#include "zebra/q_jbyt.inc"
 
#include "zebra/qstore.inc"

#ifndef CERNLIB_WINNT
      CALL UCTOH (IDHP,IDH,4,4)
#else
      CALL UCOPY (IDHP,IDH,1)
#endif
      IDN  = IDNP(1)

      L = LQT(KQT+1)
      GO TO 12

   11 L = LQ(KQS+L)
   12 IF (L.EQ.0)                  GO TO 21
      IF (IQ(KQS+L-4).NE.IDH)         GO TO 11

      IF (IDN.NE.0)  THEN
          IF (IDN.NE.IQ(KQS+L-5))     GO TO 11
        ENDIF

      LBK = L
      RETURN

C--       bank not found

   21 IF (IFLAG.EQ.0)  THEN
          LBK = 0
          RETURN
        ENDIF

#include "zebra/qtrace.inc"

      IQUEST(2)= IDH
      IQUEST(3)= IDN
      IQUEST(4)= 0
      IQUEST(5)= LQT(KQT+1)
      K = IFLAG
      IF (K.LT.100)  K=61
      CALL ZTELL (K,1)
      END
*      ==================================================
#include "zebra/qcardl.inc"
