*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZFATAL

C-    FATAL PROGRAM TERMINATION

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
C--------------    END CDE                             --------------
#if defined(CERNLIB_QMVDS)
      SAVE         INIT
#endif
#include "zebra/zfatalch.inc"

#include "zebra/zfatalre.inc"

      IF (NQERR.NE.0)              GO TO 71
      NQERR = NQERR+1
      LUN   = IQTYPE
      IF (LUN.NE.0)                GO TO 22

   21 LUN  = IQLOG
   22 IF (NQTRAC.EQ.0)             GO TO 31

C----              PRINT ZEBRA TRACE-BACK

#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_QTRHOLL)||defined(CERNLIB_A6M))
      JT = NQTRAC - 1
      WRITE (LUN,9024) MQTRAC(JT+1)
 9024 FORMAT (1X/' !!!!! ZFATAL called from ',A6)
      GO TO 28

   25 WRITE (LUN,9025) MQTRAC(JT+1)
 9025 FORMAT (14X,'called from ',A6)
   28 JT = JT - 1
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      JT = NQTRAC - 2
      WRITE (LUN,9024) MQTRAC(JT+1),MQTRAC(JT+2)
 9024 FORMAT (1X/' !!!!! ZFATAL called from ',2A4)
      GO TO 28

   25 WRITE (LUN,9025) MQTRAC(JT+1),MQTRAC(JT+2)
 9025 FORMAT (14X,'called from ',2A4)
   28 JT = JT - 2
#endif
#if defined(CERNLIB_QPRINT)
      IF (JT.GE.0)                 GO TO 25
      IF (NQFATA.EQ.0)             GO TO 49
      GO TO 41
#endif

C--                EXTERNAL CALL TO ZFATAL

   31 IF (NQFATA.NE.0)             GO TO 41
#if defined(CERNLIB_QPRINT)
      WRITE (LUN,9031)
 9031 FORMAT (1X/' !!!!! ZFATAL reached.')
#endif
      GO TO 49

C--                ZEBRA INTERNAL CALL TO ZFATAL

   41 CONTINUE
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      WRITE (LUN,9041) IQUEST(10),NQCASE
 9041 FORMAT (1X/' !!!!! ZFATAL reached from ',A6,'  for Case=',I3/1X)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      WRITE (LUN,9041) IQUEST(9),IQUEST(10),NQCASE
 9041 FORMAT (1X/' !!!!! ZFATAL reached from ',2A4,'  for Case=',I3/1X)
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_QTRHOLL))
      WRITE (LUN,9041) IQUEST(9),IQUEST(10),NQCASE
 9041 FORMAT (1X/' !!!!! ZFATAL reached from ',2A4,'  for Case=',I3/1X)
#endif
#if defined(CERNLIB_QPRINT)

      JPOS = IQBITW - 7

      DO 47  JW=11,10+NQFATA
      IT = IQUEST(JW)
      J  = JBYT (IT,JPOS,8)
      IF (J.EQ.0)                  GO TO 44
      IF (J.EQ.255)                GO TO 44

      WRITE (LUN,9043,ERR=47)  JW,IT,IT,IT
      GO TO 47

   44 WRITE (LUN,9044,ERR=47)  JW,IT,IT
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))
 9043 FORMAT (10X,'IQUEST(',I2,') = ',I9,1X,O22,1X,A6)
 9044 FORMAT (10X,'IQUEST(',I2,') = ',I9,1X,O22)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))
 9043 FORMAT (10X,'IQUEST(',I2,') = ',I9,1X,Z16,1X,A6)
 9044 FORMAT (10X,'IQUEST(',I2,') = ',I9,1X,Z16)
#endif
   47 CONTINUE

   49 WRITE (LUN,9049) JQSTOR,JQDIVI
 9049 FORMAT (1X/10X,'Current Store number =',I3,'  (JQDIVI=',I2,')')

      IF (IQVID(2).EQ.0)           GO TO 59

      WRITE (LUN,9051) IQVID
      WRITE (LUN,9052) (J,IQVREM(1,J),IQVREM(2,J),J=1,6)

 9051 FORMAT (1X/10X,'Automatic Verification Identifiers :'
     F/10X,'Current :',2X,2I11)
 9052 FORMAT (10X,'Stacked, J =',I2,' :',I6,I11,5(/22X,I2,' :',I6,I11))

   59 IF (LUN.NE.IQLOG)            GO TO 21
      NQCASE = 0
      NQFATA = 0
      CALL ZABEND

C----              RECOVERY LOOP

   71 NQERR = NQERR + 1
      IF (NQERR.GE.4)              GO TO 79
      WRITE (IQLOG,9071)
      IF (IQTYPE.EQ.0)             GO TO 79
      IF (IQTYPE.EQ.IQLOG)         GO TO 79
      WRITE (IQTYPE,9071)
 9071 FORMAT (1X/' !!!!! Stop for re-entry to ZFATAL.')
   79 CONTINUE
      CALL ABEND
      END
*      ==================================================
#include "zebra/qcardl.inc"
