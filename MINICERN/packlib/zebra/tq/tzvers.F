*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:14:56  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:28  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZVERS (IXSTOR,LBK,CHIDH,IDNP,ISELP,IFLAG)

C-    Find title bank with IDH (and IDN if non-zero)
C-    and a validity range spanning ISELP

#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
C--------------    END CDE                             -----------------
      INTEGER      CHIDH(9), IDNP(9), ISELP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HTZVE, 4HRS   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HTZVERS /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'TZVERS  ')
#endif

#include "zebra/q_jbyt.inc"

#include "zebra/qstore.inc"

      ISEL = ISELP(1)
      IDN  = IDNP(1)
      CALL UCTOH (CHIDH,IDH,4,4)

      L = LQT(KQT+1)
      GO TO 12

   11 L = LQ(KQS+L)
   12 IF (L.EQ.0)                  GO TO 21
      IF (IQ(KQS+L-4).NE.IDH)         GO TO 11

      IF (IDN.NE.0)  THEN
          IF (IDN.NE.IQ(KQS+L-5))     GO TO 11
        ENDIF

      IF (ISEL.LT.IQ(KQS+L+1))        GO TO 11
      IF (ISEL.GT.IQ(KQS+L+2))        GO TO 11

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
      IQUEST(4)= ISEL
      IQUEST(5)= LQT(KQT+1)
      K = IFLAG
      IF (K.LT.100)  K=62
      CALL ZTELL (K,1)
      END
*      ==================================================
#include "zebra/qcardl.inc"
