*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZFORM (CHID,CHFORM,IXIOP)

C-    Analyse I/O descriptor,  user called
C-    save resulting I/O characteristic tagged by CHID,
C-    return in IXIOP the index to the characteristic

#include "zebra/zvfaut.inc"
#include "zebra/zkrakc.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqform.inc"
C--------------    END CDE                             --------------
      EQUIVALENCE (NW,IQUEST(1))

      DIMENSION    IXIOP(99)
      CHARACTER    CHID*(*), CHFORM*(*)
      DIMENSION    MMID(5), MMIX(5), MMIO(5)

#if defined(CERNLIB_QMVDS)
      SAVE         MMID, MMIX, MMIO
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZFO, 4HRM   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZFORM /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZFORM  ')
#endif
#if defined(CERNLIB_QHOLL)
      DATA  MMID   / 4HQID , 2, 2, 10, 5 /
      DATA  MMIX   / 4HQIOX, 0, 0,  7, 2 /
      DATA  MMIO   / 4HQIOD, 0, 0, 50, 1 /
#endif
#if !defined(CERNLIB_QHOLL)
      DATA  MMID   / 0,      2, 2, 10, 5 /
      DATA  MMIX   / 0,      0, 0,  7, 2 /
      DATA  MMIO   / 0,      0, 0, 50, 1 /
#endif

#include "zebra/q_sbyt.inc"

#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      N = MIN (4, LEN(CHID))
      CALL UCTOH (CHID,IDH,4,N)
      IF (LQFORM.EQ.0)             GO TO 75

   12 IQCETK(121) = IDH
      LIOD  = LQ(KQSP+LQFORM-2)
      IXIOD = IQ(KQSP+LIOD+1)
      CALL MZIOCH (IQ(KQSP+LIOD+IXIOD+1),16,CHFORM)
      NW    = NW + 1
      NWIO  = IXIOD + NW
      IQ(KQSP+LIOD+1) = NWIO
      NFRIO = IQ(KQSP+LIOD-1) - NWIO

      LID = LQFORM
      IF (IDH.LT.0)  LID=LQ(KQSP+LID)
      LIX  = LQ(KQSP+LID-1)

      NWID  = IQ(KQSP+LID+1) + 1
      IQ(KQSP+LID+1)      = NWID
      IQ(KQSP+LID+NWID+3) = IDH
      IQ(KQSP+LIX+NWID)   = IXIOD
      NFRID = IQ(KQSP+LID-1) - NWID - 3

      IQUEST(2) = 64*NW + 2
      IQUEST(2) = MSBYT (IXIOD,IQUEST(2),17,15)
      IXIOP(1) = IQUEST(2)
      IF (NFRID.EQ.0)              GO TO 71
   28 IF (NFRIO.LT.16)             GO TO 73
   29 CONTINUE
#include "zebra/qtrace99.inc"
      RETURN

C------            MAKE ROOM IN THE CONTROL BANKS

C--                INCREASE ID/IX BANKS

   71 CALL MZPUSH (JQPDVS,LID,0,20,'I')
      LIX = LQ(KQSP+LID-1)
      CALL MZPUSH (JQPDVS,LIX,0,20,'I')
      GO TO 28

C--                INCREASE IOD BANK

   73 LIOD = LQ(KQSP+LQFORM-2)
      CALL MZPUSH (JQPDVS,LIOD,0,60,'I')
      GO TO 29

C--                CREATE I/O DESCRIPTOR DATA STRUCTURE

   75 CONTINUE
#if !defined(CERNLIB_QHOLL)
      CALL UCTOH ('QID ', MMID, 4,4)
      CALL UCTOH ('QIOX', MMIX, 4,4)
      CALL UCTOH ('QIOD', MMIO, 4,4)
#endif
      DO  76  J=1,2
      CALL MZLIFT (JQPDVS,L,LQFORM,1,MMID,0)
      CALL MZLIFT (JQPDVS,LIX,L,-1,MMIX,0)
   76 CONTINUE

      CALL MZLIFT (JQPDVS,L,LQFORM,-2,MMIO,0)
      IQ(KQSP+L+1) = 1
      GO TO 12
      END
*      ==================================================
#include "zebra/qcardl.inc"
