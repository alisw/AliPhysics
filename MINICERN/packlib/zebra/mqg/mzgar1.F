*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:32  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZGAR1

C-    Garbage collect division JQDIVI for not enough space
C-    System called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZGA, 4HR1   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZGAR1 /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZGAR1  ')
#endif

#include "zebra/q_sbit1.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IQVREM(1,1) = IQVID(1)
      IQVREM(2,1) = IQVID(2)
#endif

C----              Construct Memory Occupation table

      MQDVGA = 0
      MQDVWI = 0

      IF (JQDIVI.LT.3)             GO TO 24
      MQDVGA = MSBIT1 (0,JQDIVI)
      JQDVM2 = JQDIVI - JQMODE
      IF (JQDVM2.EQ.JQDVSY-1)  JQDVM2=JQDVLL
      JQDVM1 = 2
      JQSTMV = JQSTOR
      IQTNMV = 0
      IF (JQSHAR.EQ.0)             GO TO 29
      MQDVGA = MSBIT1 (MQDVGA,JQSHAR)
      GO TO 29

   24 MQDVGA = 3
      JQSTMV = -1
   29 NQDVMV = 0
      NRESAV = NQRESV
#if (defined(CERNLIB_QDEBPRI))&&(!defined(CERNLIB_QTRHOLL)||defined(CERNLIB_A6M))
      IF (NQLOGL.GE.1)  WRITE (IQLOG,9028) MQTRAC(NQTRAC-1),
     +                      JQSTOR,JQDIVI,NQRESV
 9028 FORMAT (' MZGAR1-  Auto Garbage Collection called from ',A6,
     F' for Store/Div',2I3,' Free',I7)
#endif
#if (defined(CERNLIB_QDEBPRI))&&(defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      IF (NQLOGL.GE.1)  WRITE (IQLOG,9028) MQTRAC(NQTRAC-3),
     +                      MQTRAC(NQTRAC-2),JQSTOR,JQDIVI,NQRESV
 9028 FORMAT (' MZGAR1-  Auto Garbage Collection called from ',2A4,
     F' for Store/Div',2I3,' Free',I7)
#endif
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9829) JQDIVI,JQSHAR,NQRESV
 9829 FORMAT (1X/' DEVZE MZGAR1 entered,  JQDIVI,JQSHAR,NQRESV= ',3I8)
#endif
      CALL MZTABM

C--                Construct Link Relocation table

      CALL MZTABR
      NQRESV = NQRESV + NQFREE
      IF (NQRESV.GE.0)             GO TO 51
      IF (IQPART.NE.0)             GO TO 51
      IF (JQDIVI.LT.3)             GO TO 72

C----              Shift division if not enough space

      NRESV1 = LQSTA(KQT+2) - LQEND(KQT+1) - NQMINR
      NRESV1 = MIN (NRESV1,LQEND(KQT+2)-LQ2END)

C--                Forward division

      IF (JQMODE.NE.0)             GO TO 34
      IF (JQSHAR.NE.0)  THEN
          NPOSSH = NQDMAX(KQT+JQDIVI) + NQDMAX(KQT+JQDIVN)
     +             -(LQEND(KQT+JQDIVN) - LQSTA(KQT+JQDIVI))
          GO TO 36
        ELSE
          NPOSSH = LQSTA(KQT+JQDIVI) + NQDMAX(KQT+JQDIVI)
     +             - LQSTA(KQT+JQDIVN)
          GO TO 36
        ENDIF

C--                Reverse division

   34 IF (JQSHAR.NE.0)  THEN
          NPOSSH = NQDMAX(KQT+JQDIVI) + NQDMAX(KQT+JQDIVN)
     +             -(LQEND(KQT+JQDIVI) - LQSTA(KQT+JQDIVN))
        ELSE
          NPOSSH = LQEND(KQT+JQDIVN)
     +             - (LQEND(KQT+JQDIVI) - NQDMAX(KQT+JQDIVI))
        ENDIF

   36 NSH = (LQEND(KQT+JQDIVI)-LQSTA(KQT+JQDIVI)) / 8
      NSH = MAX (NSH,24) - NQRESV
      NSH = MIN (NSH, NPOSSH, NRESV1)

      IF (NSH+NQRESV.LT.0)         GO TO 72

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.7)
     +WRITE (IQLOG,9836) JQDIVI,JQSHAR,NQRESV
     +,           NRESV1,NPOSSH,NSH
     +,           JQGAPM,JQGAPR
 9836 FORMAT (' DEVZE MZGAR1,  JQDIVI,JQSHAR,NQRESV=    ',3I8/
     F16X,'NRESV1,NPOSSH,NSH=',3I8/
     F16X,'JQGAPM,JQGAPR=    ',2I8)
#endif
      NQRESV = NQRESV + NSH
      NQDVMV = - NSH
      CALL MZTABS

C----              Relocate + memory move

   51 NWIN = NQRESV - NRESAV
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.1)  WRITE (IQLOG,9051) NWIN,NQDVMV
 9051 FORMAT (10X,'Wins',I7,' words, Shift by',I7)
#endif
      CALL MZTABX
      CALL MZTABF
      IF   (NQNOOP)          68, 53, 67
   53 CALL MZGSTA (NQDGAF(KQT+1))
      CALL MZRELX
   67 CALL MZMOVE
   68 IF (NQRESV.LT.0)             GO TO 71
#include "zebra/qtrace99.inc"
      RETURN

C--------          Not enough space

   71 IF (IQPART.NE.0)             GO TO 29
   72 IQUEST(11) = NQRESV
      IQUEST(12) = JQSTOR
      IQUEST(13) = JQDIVI
#if defined(CERNLIB_QDEBPRI)
      IF (NQLOGL.GE.1)  WRITE (IQLOG,9072) NQRESV
 9072 FORMAT (10X,'Not enough space, Free',I7)
#endif
      IF (NQPERM.NE.0)             GO TO 999
      IF (JQKIND.NE.1)             GO TO 91
      CALL ZTELL (99,1)

C------            Error conditions

   91 NQCASE = 1
      NQFATA = 1
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
