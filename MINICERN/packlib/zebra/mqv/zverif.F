*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVERIF (IXSTP,IFRETN,TEXTID)

C-    Verify integrity of store IXVSTO
C-        IXVSTO = -1 : verify all stores
C-                 -2 : verify stores selected in the ZVAUTO mask

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/zvfaut.inc"
#include "zebra/zvfwkc.inc"
C--------------    END CDE                             --------------
      CHARACTER    TEXTID*(*)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZVER, 4HFY   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZVERIF /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZVERIF  ')
#endif

#include "zebra/q_jbit.inc"

#include "zebra/qtraceq.inc"

      TEXTHD = TEXTID
      CHWARN = ' !!!!     '
      CALL VZERO (IFLHD,8)
      JVMODE = -1

C--                Check the process parameters

      CALL ZVDO1

C----              Select the store to be verified

      IXVSTO = IXSTP
      IF (IXVSTO.LT.0)             GO TO 12
      NLOOP  = 0
      CALL MZSDIV (IXVSTO,0)
      JDVINI = JQDIVI
      JST    = JQSTOR
      GO TO 19

   12 NLOOP  = 7
      JDVINI = 0
      JST    = 0
      IF (IXVSTO.EQ.-1)            GO TO 19

      IF (IXVSTO.NE.-2)    CALL ZFATAM ('not a valid store index.')

      NLOOP = -7
      JST   = -1
      GO TO 17

   16 IF (NLOOP.EQ.0)              GO TO 71
      JST = JQSTOR
   17 JST = JST + 1
      IF (JST.GT.NQSTOR)           GO TO 71
      IF (NLOOP.GE.0)              GO TO 19
      IF (JBIT(MASKST,JST+1).EQ.0) GO TO 17

   19 IXVSTO = 0
      CALL SBYT (JST,IXVSTO,27,4)
      IF (NQALLO(JST+1).LT.0)      GO TO 16

C------            Do the next store

      CALL MZSDIV (IXVSTO,-1)

      CALL VZERO (JVMODE,7)
      IF (NQLOGL.GE.2)  CALL ZVFPRI

C--                Check the store parameters

      CALL ZVDO2

      LSTOLO = LQSTA(KQT+1)
      LSTOHI = LQEND(KQT+20)

C--                If check of single division

      JQDIVI = JDVINI
      IF (JDVINI.NE.0)             GO TO 22

C----              Check all link areas

      JVMODE = 1
      L = LQSYSS(KQT+1)
      IF (L.NE.0)  THEN
          IQ(KQS+L+3) = IQ(KQS+L+2) + NQLINK
          CALL ZVDOLA (IQ(KQS+L+1))
        ENDIF

C------            Check all divisions

      JQDIVI = 0
   21 IF (JQDIVI.EQ.20)            GO TO 16
      IF (JDVINI.NE.0)             GO TO 71
      IF (JQDIVI.EQ.JQDVLL)  THEN
          JQDIVI = JQDVSY
        ELSE
          JQDIVI = JQDIVI + 1
        ENDIF

C----              Do next division

   22 JVMODE = 2
      LDIVLO = LQSTA(KQT+JQDIVI)
      LDIVHI = LQEND(KQT+JQDIVI)
      MASKTO = IQRCU(KQT+JQDIVI)
      IF (LDIVLO.GE.LDIVHI)        GO TO 21

      NAMEPR(1) = IQDN1(KQT+JQDIVI)
      NAMEPR(2) = IQDN2(KQT+JQDIVI)
      IFLDV = 0
      LBKNX = LDIVLO
      IFLBK = 7

C--                Next bank

   24 IF (LBKNX.GE.LDIVHI)         GO TO 21

      CALL MZCHLN (-7,LBKNX)
      IF (IQFOUL.NE.0)             GO TO 41

      LBKNX = IQNX
      IF (IQND.LT.0)               GO TO 24

      LBKLN = IQLN
      LBKLS = IQLS
      IBKXX = JBIT(IQ(KQS+LBKLS),IQDROP)
      IDBK(1)  = IQ(KQS+IQLS-4)
      IDBK(2)  = IQ(KQS+IQLS-5)
      NBKNN(1) = IQNIO
      NBKNN(2) = IQNL
      NBKNN(3) = IQNS
      NBKNN(4) = IQND

      IFLBK = 0
      CALL ZVDOBK
      IF (LBKNX.LE.LDIVHI)         GO TO 24

      NFATAL = NFATAL + 1
      CALL ZVFPRI
      N = LDIVHI - LBKNX
      WRITE (IQLOG,9027) CHWARN,N
 9027 FORMAT (A,'Bank overshoots the division end by',I7,' words')
      GO TO 21

C----     bank chaining clobbered

   41 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9041) CHWARN,LBKNX
 9041 FORMAT (A,'Bank chaining clobbered at adr',I10)

      L = LBKNX
   44 L = L + 1
      IF (L.GE.LDIVHI)             GO TO 47

      CALL MZCHLN (-7,L)
      IF (IQFOUL.NE.0)             GO TO 44
      LR = L

   46 CALL MZCHLN (-7,IQNX)
      IF (IQFOUL.NE.0)             GO TO 44
      IF (IQND.LT.0)               GO TO 46
      L = LR

   47 WRITE (IQLOG,9047) L
 9047 FORMAT (10X,'recover at adr',I10)

      LBKNX = L
      GO TO 24

C------            Finished

   71 IQUEST(1) = NFATAL
      IF (NFATAL+NWARN.EQ.0)       GO TO 999

      WRITE (IQLOG,9071) CHWARN, NFATAL,NWARN
 9071 FORMAT (A,'ZVERIF found',I5,' fatal and',I5,' warning conditions')

      IF (NFATAL.EQ.0)             GO TO 999
      IF (JDVINI.NE.0)  THEN
          IF (IFRETN.NE.0)         GO TO 999
        ENDIF

      CALL SBYT (LFAILS,IXVSTO,27,4)
      CALL MZSDIV (IXVSTO,-1)
      JQDIVI = LFAILD
      IQVSTA = 0
      CALL ZFATAM ('trouble in ZVERIF.')
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
