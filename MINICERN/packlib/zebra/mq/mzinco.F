*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZINCO (LIST)

C-    Initialise all passive COMMONs, normally called from MZEBRA
C-        but it may be user called for non-ZEBRA applications

#include "zebra/zbcd.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zceta.inc"
#include "zebra/zheadp.inc"
#include "zebra/zmach.inc"
#include "zebra/znatur.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/quest.inc"
C--------------    End CDE                             --------------
      DIMENSION    LIST(9)

#include "zebra/q_jbit.inc"


C--                Clear /ZSTATE/

      CALL VZERO (IQUEST,100)
      CALL VZERO (IQVID,18)
      CALL VZERO (NQPHAS,15)
#include "zebra/qversion.inc"

C----              Ready /ZMACH/

      NQBITW = IQBITW
      NQBITC = IQBITC
      NQCHAW = IQCHAW
      NQLNOR = 58
      NQLMAX = 58
      NQLPTH =  0
      NQRMAX = 132
      IQLPCT = IQBLAN
      IQNIL  = 16744448
#if defined(CERNLIB_CRAY)
      IQNIL  = 0777770516040020000000B
#elif defined(CERNLIB_CDC)
      IQNIL  = O"17770516040000200000"
#endif


C----              Ready /ZBCD/ and /ZBCDCH/

      CQALLC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-*/()$= ,.'
      CQALLC(65:90) = 'abcdefghijklmnopqrstuvwxyz'
      CQALLC(48:64) = '#''!:"_]&@?[>< ^;%'
      CQALLC(91:96) = '{|}~`?'
#if defined(CERNLIB_CDC)
      CQALLC(91:96) = '??????'
#endif
#if defined(CERNLIB_QEBCDIC)
      CQALLC(61:61) = CHAR(224)
#endif
#if !defined(CERNLIB_QEBCDIC)
      CQALLC(61:61) = CHAR(92)
#endif
      CALL UCTOH1 (CQALLC, IQLETT, 96)
      CALL UCTOH1 (' 1234567890', IQNUM2, 11)

      CALL IZHNUM (IQLETT,NQHOLL,95)
      NQHOL0 = NQHOLL(45)

C----              READY  /ZCETA/
C--   Table entry IQCETA(JH+1) contains the CETA value for
C-    the character of internal representation JH

      CALL VFILL (IQCETA,NQTCET,96)

      DO 24  JC=95,1,-1
      JH = NQHOLL(JC)
   24 IQCETA(JH+1) = JC

#if (defined(CERNLIB_QEBCDIC))&&(!defined(CERNLIB_CERNWYL))
      IQCETA(1+ 64) = 45
      IQCETA(1+189) = 54
      IQCETA(1+173) = 58
      IQCETA(1+224) = 61
      IQCETA(1+139) = 91
      IQCETA(1+192) = 91
      IQCETA(1+155) = 93
      IQCETA(1+208) = 93
#endif
#if (defined(CERNLIB_QEBCDIC))&&(defined(CERNLIB_CERNWYL))
      IQCETA(1+ 64) = 45
      IQCETA(1+189) = 54
      IQCETA(1+173) = 58
      IQCETA(1+224) = 61
      IQCETA(1+139) = 91
      IQCETA(1+192) = 91
      IQCETA(1+155) = 93
      IQCETA(1+208) = 93
      IQCETA(1+ 95) = 94
      IQCETA(1+161) = 94
#endif
#if defined(CERNLIB_QCDCODE)
      IQCETA(1) = 51
#endif
C--   Table IQTCET(JH+1) is like IQCETA but for 6-bit packing

      DO 26  JL=1,NQTCET
      J = IQCETA(JL)
      IF (J.GE.64)  THEN
        IF (J.GE.94)  THEN
          J = 57
        ELSEIF (J.EQ.93)  THEN
          J = 42
        ELSEIF (J.EQ.92)  THEN
          J = 40
        ELSEIF (J.EQ.91)  THEN
          J = 41
        ELSEIF (J.EQ.64)  THEN
          J = 51
        ELSE
C--       lower case mapped to upper case
          J = J - 64
        ENDIF
       ENDIF
   26 IQTCET(JL) = J

C----              Ready /ZNATUR/

      QPI    = 4.*ATAN(1.)
      QPI2   = 2.*QPI
      QPIBY2 = QPI/2.
      QPBYHR = .0002998

C----              Ready COMMON  /ZUNIT/

#include "mzeunit.inc"
      IQLOG  = IQPRNT
#include "mzeunit2.inc"
      ITYPE = IQTYPE
      IF (ITYPE.EQ.0)   ITYPE = IQLOG
      NLIST  = LIST(1)
      IF    (NLIST)          32, 38, 33
   32 NLIST = -NLIST
      IF (JBIT(NLIST,2).NE.0)  NQLOGD = -2
      IF (JBIT(NLIST,1).NE.0)  IQLOG = ITYPE
      IQPRNT = IQLOG
      GO TO 38

   33 NQLOGD = LIST(2)

      IF (NLIST.EQ.1)              GO TO 38
      IF     (LIST(3).NE.0)   THEN
          IF (LIST(3).LT.0)       THEN
              IQLOG = ITYPE
            ELSE
              IQLOG = LIST(3)
            ENDIF
        ENDIF
      IQPRNT = IQLOG

      IF (NLIST.EQ.2)              GO TO 38
      IF     (LIST(4).NE.0)   THEN
          IF (LIST(4).LT.0)        THEN
              IQPRNT = ITYPE
            ELSE
              IQPRNT = LIST(4)
            ENDIF
        ENDIF

   38 IQPR2  = IQPRNT
      NQLOGM = NQLOGD


      IQDLUN = IQPRNT
      IQFLUN = IQPRNT
      IQHLUN = IQPRNT
      NQUSED = 0

C----              Ready COMMON  /ZHEADP/

      CALL VBLANK (IQHEAD,20)
      CALL VZERO  (IQDATE,7)
      CALL DATIME (IQDATE,IQTIME)

#if defined(CERNLIB_CDC)
      CALL XSETIO
#endif

      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
