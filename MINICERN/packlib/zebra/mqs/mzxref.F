*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:12  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZXREF (IXFRP,IXTOP,CHOPT)

C-    Set cross-reference division IXFR pointing to IXTO

#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      DIMENSION    IXFRP(9), IXTOP(9)
      CHARACTER    *(*)  CHOPT
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZXR, 4HEF   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZXREF /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZXREF  ')
#endif

#include "zebra/q_jbyt.inc"
#include "zebra/q_sbit0.inc"
#include "zebra/q_sbit1.inc"


      IXFR = IXFRP(1)
      IXTO = IXTOP(1)

#include "zebra/qtrace.inc"

      CALL UOPTC (CHOPT,'ARC',IQUEST)
      MODE = IQUEST(1)
      IF (IQUEST(2).NE.0)  MODE=-1
      IF (IQUEST(3).NE.0)  MODE=-2

C-        MODE  +1 Add, 0 set, -1 Remove, -2 Contained

      CALL MZSDIV (IXFR,0)
      IF (JQDIVI.EQ.0)             GO TO 91
      IF (MODE.EQ.-2)              GO TO 27

      JST = JBYT (IXTO,27,6)
      IF (JST-16.EQ.JQSTOR)        GO TO 31
      IF    (JST.NE.JQSTOR)        GO TO 92

C--                SINGLE DIVISION INDEX

      JDV = JBYT (IXTO,1,26)
      IF (JDV.GE.25)               GO TO 93
      IF (JDV.EQ.0)                GO TO 93
      IQUEST(1) = JDV
      NBI = 1
      GO TO 34

C--                CONTAINED DIVISION

   27 IQRCU(KQT+JQDIVI) = 0
      IQRTO(KQT+JQDIVI) = 0
      IQRNO(KQT+JQDIVI) = 0
      GO TO 999

C--                COMPOSITE INDEX

   31 CALL UBITS (IXTO,26,IQUEST,NBI)

      IF (NBI.EQ.0)                GO TO 93
      IF (IQUEST(NBI).GE.25)       GO TO 93

   34 IF (MODE.NE.0)               GO TO 41
      IQRTO(KQT+JQDIVI) = 0
      IQRNO(KQT+JQDIVI) = 9437183

   41 DO 49  JBI=1,NBI
      JDV = IQUEST(JBI)
      IF (MODE.LT.0)               GO TO 47

C--       add
      IQRTO(KQT+JQDIVI) = MSBIT1 (IQRTO(KQT+JQDIVI),JDV)
      IF (JDV.GE.21 .AND. JDV.LT.24)   GO TO 49
      IQRNO(KQT+JQDIVI) = MSBIT1 (IQRNO(KQT+JQDIVI),JDV)
      GO TO 49

C--       remove
   47 IQRTO(KQT+JQDIVI) = MSBIT0 (IQRTO(KQT+JQDIVI),JDV)
      IQRNO(KQT+JQDIVI) = MSBIT0 (IQRNO(KQT+JQDIVI),JDV)
   49 CONTINUE
      CALL MZXRUP
#include "zebra/qtrace99.inc"
      RETURN

C------            ERROR CONDITIONS

   93 NQCASE = 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = 3
      IQUEST(11) = IXFR
      IQUEST(12) = IXTO
      IQUEST(13) = MODE
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
