*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:43  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZSHUNT (IXSTOR,LSHP,LSUPP,JBIASP,IFLAGP)

C-    RE-CONNECT BANK OR LINEAR D/S, USER CALLED

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
C--------------    END CDE                             --------------
      DIMENSION    LSHP(9),LSUPP(9),JBIASP(9),IFLAGP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HZSHU, 4HNT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HZSHUNT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'ZSHUNT  ')
#endif

#include "zebra/q_jbyt.inc"
#include "zebra/q_locf.inc"

#include "zebra/qtraceq.inc"

      LSH = LSHP(1)
      IF (LSH.EQ.0)                GO TO 999
      LSUP  = LSUPP(1)
      JBIAS = JBIASP(1)
      IFLAG = IFLAGP(1)

#include "zebra/qstore.inc"

#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LSH)
      IF (IQFOUL.NE.0)             GO TO 91
#endif
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.2)  THEN
          IF (JBIAS.GE.2)  LSUP=0
          WRITE (IQLOG,9011) JQSTOR,LSH,LSUP,JBIAS,IFLAG,IQID
        ENDIF
 9011 FORMAT (' ZSHUNT-  Store',I3,' LSH/LSUP/JBIAS/IFLAG='
     F,2I9,1X,I6,1X,I3,' IDH= ',A4)
#endif

C----              LINKS AT EXIT POINT

      KEX  = LQ(KQS+LSH+2)
      LNEX = LQ(KQS+LSH)

C----              LINKS AT INSERTION POINT

      LPRE = 0
      IF   (JBIAS-1)         21, 25, 28
   21 CONTINUE
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LSUP)
      IF (IQFOUL.NE.0)             GO TO 92
      IF (IQNS+JBIAS.LT.0)         GO TO 93
#endif
      KIN  = LSUP + JBIAS
      LNIN = LQ(KQS+KIN)
      LUP  = LSUP
      IF (JBIAS.NE.0)              GO TO 29
      LPRE = LUP
      LUP  = LQ(KQS+LUP+1)
      GO TO 29

   25 LNIN = LSUP
      IF (LNIN.EQ.0)               GO TO 26
#if defined(CERNLIB_QDEBUG)
      CALL MZCHLS (-7,LSUP)
      IF (IQFOUL.NE.0)             GO TO 92
#endif
      KIN  = LQ(KQS+LNIN+2)
      LUP  = LQ(KQS+LNIN+1)
      GO TO 29

   26 KIN  = LOCF(LSUPP(1)) - LQSTOR
      LUP  = 0
      GO TO 29

   28 KIN  = 0
      LNIN = 0
      LUP  = 0
      IF (KEX.EQ.0)                GO TO 51

C--                CHECK INSERT POINT = EXIT POINT

   29 IF (KIN.EQ.KEX)              GO TO 999

C--                CHECK LINEAR STRUCTURE CONTAINED

#if defined(CERNLIB_QDEBUG)
      L = MAX  (LNIN,LPRE)
      IF (L.EQ.0)                  GO TO 51

      IF (L.GE.LQEND(KQT+20))         GO TO 94
      IF (L.GE.LQEND(KQT+JQDVLL))     GO TO 43
      JQDIVI = 2
      IF (L.GE.LQEND(KQT+2))          GO TO 44
      IF (L.GE.LQSTA(KQT+2))          GO TO 45
      JQDIVI = 1
      GO TO 45

   43 JQDIVI = JQDVSY - 1
   44 JQDIVI = JQDIVI + 1
      IF (L.GE.LQEND(KQT+JQDIVI))     GO TO 44

   45 IF (LSH.LT.LQSTA(KQT+JQDIVI))   GO TO 94
      IF (LSH.GE.LQEND(KQT+JQDIVI))   GO TO 94
#endif

C----              SHUNT LINEAR STRUCTURE

   51 IF (LNEX.EQ.0)               GO TO 58
      IF (IFLAG.EQ.0)              GO TO 57
#if defined(CERNLIB_QDEBUG)
      L  = LSH
   53 CALL MZCHLS (-7,LNEX)
      IF (IQFOUL.NE.0)             GO TO 95
      L  = LNEX
      LNEX = LQ(KQS+LNEX)
      IF (LNEX.NE.0)               GO TO 53
#endif
      LNEX = LSH
   55 LEND = LNEX
      LQ(KQS+LEND+1) = LUP
      LNEX = LQ(KQS+LEND)
      IF (LNEX.NE.0)               GO TO 55
      GO TO 71

C----              SHUNT SINGLE BANK

   57 CONTINUE
#if defined(CERNLIB_QDEBUG)
      L  = LSH
      CALL MZCHLS (-7,LNEX)
      IF (IQFOUL.NE.0)             GO TO 95
#endif
   58 LEND = LSH
      LQ(KQS+LSH+1) = LUP

C----              CONNECTIONS

C--                BRIDGE OLD POSITION

   71 IF (KEX .NE.0)  LQ(KQS+KEX) = LNEX
      IF (LNEX.NE.0)  LQ(KQS+LNEX+2) = KEX

C--                CONNECT START

      IF (KIN.NE.0)  THEN
          LQ(KQS+KIN) = LSH
        ELSE
          LSUPP(1) = LSH
        ENDIF
      LQ(KQS+LSH+2) = KIN

C--                CONNECT END

      LQ(KQS+LEND) = LNIN
      IF (LNIN.NE.0)  LQ(KQS+LNIN+2) = LEND

#include "zebra/qtrace99.inc"
      RETURN

C------            ERROR CONDITIONS

#if defined(CERNLIB_QDEBUG)
   95 NQCASE = 1
      NQFATA = 1
      IQUEST(16) = LNEX
   94 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 1
      IQUEST(15) = L
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
#endif
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 4
      IQUEST(11) = LSH
      IQUEST(12) = LSUP
      IQUEST(13) = JBIAS
      IQUEST(14) = IFLAG
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
