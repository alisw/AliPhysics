*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:24  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVDO2

C-    Subsidiary to ZVERIF, check store table

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/zvfwkc.inc"
C--------------    End CDE                             --------------
      PARAMETER   (NVECPR=342)
      DIMENSION    MVECPR(NVECPR)
      EQUIVALENCE (MVECPR,NQSTOR)
      CHARACTER    TEXT*8
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZVDO, 4H2    /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZVDO2  /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZVDO2   ')
#endif


      NQCASE = 0
      JDV    = 0
      IF (LQSTOR.NE.LQASTO+KQS)    GO TO 81
      IF (NQFEND.LE.0)             GO TO 82

C--       check monotony of store parameters

      IF (NQSTRU.GT.NQREF)         GO TO 83
      IF (NQREF.GT.NQLINK)         GO TO 84
      LAST = NQLINK

      JDV = 1
   22 IF (JDV.EQ.JQDVLL+1)  JDV = JQDVSY

      IF (LQSTA(KQT+JDV).LT.LAST)     GO TO 85
      LAST = LQSTA(KQT+JDV)
      IF (LQEND(KQT+JDV).LT.LAST)     GO TO 86
      LAST = LQEND(KQT+JDV)
      JDV  = JDV + 1
      IF (JDV.LE.20)               GO TO 22

C--       check the fence is intact

      DO 24  J=1,NQFEND
      IF (LQ(KQS+J-NQFEND).NE.IQNIL) GO TO 61
   24 CONTINUE
      RETURN

C----              Fence overwritten

   61 WRITE (IQLOG,9061)
      TEXT = 'FENCE'
      WRITE (IQLOG,9099) TEXT, (LQ(KQS+J-NQFEND),J=1,NQFEND)
      GO TO 98

C--                Non-monotony

   86 NQCASE = 1
   85 NQCASE = NQCASE + 1
   84 NQCASE = NQCASE + 1
   83 NQCASE = NQCASE + 1
   82 NQCASE = NQCASE + 1
   81 NQCASE = NQCASE + 1
      NQFATA = 1
      IQUEST(1) = JDV

      WRITE (IQLOG,9091)
      TEXT = 'IQTABV'
      WRITE (IQLOG,9099) TEXT,(IQTABV(KQT+J),J=1,16)
      TEXT = 'LQSYSS'
      WRITE (IQLOG,9099) TEXT,(LQSYSS(KQT+J),J=1,10)
      TEXT = 'LQSYSR'
      WRITE (IQLOG,9099) TEXT,(LQSYSR(KQT+J),J=1,10)
      TEXT = 'IQTDUM'
      WRITE (IQLOG,9099) TEXT,(IQTDUM(KQT+J),J=1,22)
      TEXT = 'LQSTA'
      WRITE (IQLOG,9099) TEXT,(LQSTA(KQT+J),J=1,21)
      TEXT = 'LQEND'
      WRITE (IQLOG,9099) TEXT,(LQEND(KQT+J),J=1,20)
      TEXT = 'NQDMAX'
      WRITE (IQLOG,9099) TEXT,(NQDMAX(KQT+J),J=1,20)
      TEXT = 'IQMODE'
      WRITE (IQLOG,9099) TEXT,(IQMODE(KQT+J),J=1,20)
      TEXT = 'IQKIND'
      WRITE (IQLOG,9099) TEXT,(IQKIND(KQT+J),J=1,20)
      TEXT = 'IQRCU'
      WRITE (IQLOG,9099) TEXT,(IQRCU(KQT+J),J=1,20)
      TEXT = 'IQRTO'
      WRITE (IQLOG,9099) TEXT,(IQRTO(KQT+J),J=1,20)
      TEXT = 'IQRNO'
      WRITE (IQLOG,9099) TEXT,(IQRNO(KQT+J),J=1,20)
      TEXT = 'NQDINI'
      WRITE (IQLOG,9099) TEXT,(NQDINI(KQT+J),J=1,20)
      TEXT = 'NQDSIZ'
      WRITE (IQLOG,9099) TEXT,(NQDSIZ(KQT+J),J=1,20)

   98 CONTINUE
#include "zebra/qtraceq.inc"
#include "zebra/qtofatal.inc"

 9061 FORMAT (1X/' ZVDO2.   Fence overwritten, dump follows'/1X)
 9091 FORMAT (1X/' ZVDO2.   Store parameter table overwritten'
     F/10X,'dump the /MZCC/ equivalent for the current store'/1X)
#if defined(CERNLIB_B32)
 9099 FORMAT (3X,A,1X,5Z10/(12X,5Z10))
#endif
#if !defined(CERNLIB_B32)
 9099 FORMAT (3X,A,1X,4Z18/(12X,4Z18))
#endif
      END
*      ==================================================
#include "zebra/qcardl.inc"
