*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:35  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZINQD (IXDIV)

C-    Return the attributes of division IXDIV :
C-        IQUEST(1)  mode = 0 forward, 1 reverse division
C-        IQUEST(2)  kind = 1 short range, 2 long, 4 package, 8 system
C-        IQUEST(3)  start adr LSTA
C-        IQUEST(4)  end   adr LEND
C-        IQUEST(5)  if mode=0 : LEMAX highest possible value of LEND
C-                   if mode=1 : LSMIN  lowest possible value of LSTA
C-        IQUEST(6)  initial size
C-        IQUEST(7)  declared maximum size
C-        IQUEST(8)  char 1:4 of BCD name of division
C-        IQUEST(9)  char 5:8
C-        IQUEST(10) JDVSH, if non-zero: sharing neighbour division

C-    ie.:   memory currently occupied =  LEND - LSTA
C-
C-           immediate free space  = LEMAX - LEND  if MODE = 0
C-                                 = LSTA - LSMIN  if MODE = 1

#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZIN, 4HQD   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZINQD /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZINQD  ')
#endif

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      CALL MZSDIV (IXDIV,4)

      MODE = JBIT (IQMODE(KQT+JQDIVI), 1)
      KIND = JBYT (IQKIND(KQT+JQDIVI),21,4)
      LSTA = LQSTA(KQT+JQDIVI)
      LEND = LQEND(KQT+JQDIVI)
      NMAX = NQDMAX(KQT+JQDIVI)
      IF (JQDIVI.LT.3)             GO TO 41
      IF (MODE.NE.0)               GO TO 31

C--                FORWARD DIVISION

      JDVNE = JQDIVI + 1
      IF (JQDIVI.EQ.JQDVLL)        GO TO 24
      IF (JQDIVI.EQ.20)            GO TO 25

      IF (JBYT(IQKIND(KQT+JDVNE),21,4).NE.KIND)  GO TO 25
      IF (JBIT(IQMODE(KQT+JDVNE),1)   .EQ.MODE)  GO TO 25
      JDVSH = JDVNE
      LMAX  = LQSTA(KQT+JDVNE)
      GO TO 78

   24 JDVNE = JQDVSY
   25 LMAX  = MIN (LSTA+NMAX, LQSTA(KQT+JDVNE) )
      JDVSH = 0
      GO TO 78

C--                REVERSE DIVISION

   31 JDVNE = JQDIVI - 1
      IF (JQDIVI.EQ.JQDVSY)        GO TO 34
      IF (JBYT(IQKIND(KQT+JDVNE),21,4).NE.KIND)  GO TO 35
      IF (JBIT(IQMODE(KQT+JDVNE),1)   .EQ.MODE)  GO TO 35
      JDVSH = JDVNE
      LMAX  = LQEND(KQT+JDVNE)
      GO TO 78

   34 JDVNE = JQDVLL
   35 LMAX  = MAX (LEND-NMAX, LQEND(KQT+JDVNE) )
      JDVSH = 0
      GO TO 78

C----              DIVISIONS 2 + 1

   41 IF (JQDIVI.EQ.1)             GO TO 44
      LMAX  = LQEND(KQT+1) + NQMINR
      JDVSH = 1
      GO TO 78

   44 LMAX  = LQSTA(KQT+2) - NQMINR
      JDVSH = 2

C----              Store vector

   78 IQUEST(1)  = MODE
      IQUEST(2)  = KIND
      IQUEST(3)  = LSTA
      IQUEST(4)  = LEND
      IQUEST(5)  = LMAX
      IQUEST(6)  = NQDINI(KQT+JQDIVI)
      IQUEST(7)  = NMAX
      IQUEST(8)  = IQDN1(KQT+JQDIVI)
      IQUEST(9)  = IQDN2(KQT+JQDIVI)
      IQUEST(10) = JDVSH
      IQUEST(11) = 0
      IQUEST(12) = 0
#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.2)
     +     WRITE (IQLOG,9071) JQDIVI,JQSTOR,(IQUEST(J),J=1,10)

 9071 FORMAT (/' MZINQD.  Inquiry for division',I3,' in store',I3,
     F' returns :'
     F/10X,'mode kind     LSTA     LEND    LEXTR  inisize  maxsize'
     F,'  name1/2 JDVSH'
     F/9X,2I5,5I9, 2X,2A4,I5)
#endif
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
