*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZXRPR (IXSTOR)

C-    Print cross-reference table for store IXSTOR

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
      CHARACTER    CHWK*24, CHYN(6)*1
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZXR, 4HPR   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZXRPR /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZXRPR  ')
#endif
      DATA  CHYN   / 'n', 'y', '.', 's', 'u', '.' /

#include "zebra/q_jbit.inc"

#include "zebra/qtrace.inc"

      JSTINI = JQSTOR
      JDVINI = JQDIVI
      IF (IXSTOR.EQ.-7)  THEN
          IF (JQDIVI.EQ.0)         RETURN
          JDV  = JQDIVI
          IFLS = 7
          GO TO 22
        ENDIF

      IFLS = 0
      CALL MZSDIV (IXSTOR,0)
      MASKDV = MZDVAC (IXSTOR)
      IF (MASKDV.EQ.0)  MASKDV = 1048575

      WRITE (IQLOG,9001) JQSTOR
 9001 FORMAT (1X/' MZXRPR.  Cross reference table for store',I3)

C----              Loop all divisions, if not specific

      JDV = 1
   21 IF (JDV.EQ.JQDVLL+1)  JDV = JQDVSY
      IF (JBIT(MASKDV,JDV).EQ.0)   GO TO 49

      WRITE (IQLOG,9007)
 9007 FORMAT (1X)

C--                Do next division

   22 WRITE (IQLOG,9022) JDV
 9022 FORMAT (10X,'Table entries for division',I3)

      CALL UPKBYT (IQRCU(KQT+JDV),1,IQUEST(1),24,0)
      DO 42  J=1,24
      I = IQUEST(25-J)
   42 CHWK(J:J) = CHYN(I+1)
      WRITE (IQLOG,9042) CHWK(5:9),
     +                   CHWK(10:14),CHWK(15:19),CHWK(20:24)

      CALL UPKBYT (IQRTO(KQT+JDV),1,IQUEST(1),24,0)
      DO 43  J=1,24
      I = IQUEST(25-J)
   43 CHWK(J:J) = CHYN(I+3)
      WRITE (IQLOG,9043) CHWK(1:4),CHWK(5:9),
     +                   CHWK(10:14),CHWK(15:19),CHWK(20:24)

      CALL UPKBYT (IQRNO(KQT+JDV),1,IQUEST(1),24,0)
      DO 44  J=1,24
      I = IQUEST(25-J)
   44 CHWK(J:J) = CHYN(I+5)
      WRITE (IQLOG,9044) CHWK(5:9),
     +                   CHWK(10:14),CHWK(15:19),CHWK(20:24)

 9042 FORMAT (12X,'active: ',5X,4(1X,A))
 9043 FORMAT (12X,'   set: ',   5(1X,A))
 9044 FORMAT (12X,' unset: ',5X,4(1X,A))

      IF (IFLS.NE.0)               GO TO 999
   49 JDV = JDV + 1
      IF (JDV.LE.20)               GO TO 21
      IF (JQSTOR.EQ.JSTINI)  JQDIVI = JDVINI

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
