*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZIDIA

C-    Diagnostic handling routine for FZIN and FZCOPY

C-    The error return code JRETCD is now set :

C-    -4  all segments to be skipped : FZIMTB
C-    -3  trying to read beyond EoD
C-    -2  normal start/end of run : FZIFFN, FZIFFX
C-    -1  normal EoF : FZIFFN, FZIREC
C-     1  unexpected EoF :  FZIFFN, FZIPHR, FZIPHA
C-     2  unexpected start/end of run : FZIFFN
C-     3  not enough memory : FZIMTB
C-     4  Bad User Handling : FZIN, FZIFFX, FZIMTB
C-     5  Bad Data : FZIFFN, FZIFFX, FZIPHR, FZIPHA, FZIASC, FZIREL
C-     6  Bad Construction : FZIFFN, .FFX, .REC, .PHR, .PHM, .PHA
C-     7  Read Error : FZIFFN, FZIPHR, FZIASC
C-     8  Emergency Stop

C-    The old scheme was :

C-     1  unexpected EoF :  FZIFFN, FZIPHR, FZIPHA
C- -4     also : all segments to be skipped in FZIMTB
C-  5            bad data in FZIREL
C- -1  2  normal EoF : FZIFFN, FZIREC
C-  4     also : user error handling the segment table in FZIMTB
C-  2  3  unexpected start/end of run : FZIFFN
C-  3     also : not enough memory in FZIMTB
C- -2  4  normal start/end of run : FZIFFN, FZIFFX
C-     5  Bad Data : FZIFFN, FZIFFX, FZIPHR, FZIPHA, FZIASC
C-     6  Bad Construction : FZIFFN, .FFX, .REC, .PHR, .PHM, .PHA
C-     7  Read Error : FZIFFN, FZIPHR, FZIASC
C-  4  8  user routine / memory not connected : FZIFFX
C-  8  -  Emergency Stop

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (LRTYP,IDI(2))
      CHARACTER    ERRMSG*16
      CHARACTER    TEXT*4
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZID, 4HIA   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZIDIA /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZIDIA  ')
#endif

#include "zebra/q_jbit.inc"

#include "zebra/qtrace.inc"

      IF (JRETCD.GE.0)             GO TO 40

C-----------------------------------------------------
C-                 Normal end-of-file / end-of-run
C-----------------------------------------------------

      J = -JRETCD
      GO TO ( 24, 31, 21), J

C--                Reading beyond end-of-data

   21 CONTINUE
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.-2)  WRITE (IQLOG,9021) LUNI
 9021 FORMAT (1X/' FZIDIA.  LUN=',I4,'  Attempt to read beyond E-o-D')
#endif
      IF (IACTVI.EQ.7)       CALL ZFATAM ('FZIDIA - reading beyond EOD')
      IACTVI = 7
      GO TO 27

C--                Normal EOF

   24 N = IQ(KQSP+LQFI+11) + 1
      IQ(KQSP+LQFI+11) = N
      IACTVI = MAX (IACTVI,4) + 1
      IF (JBIT(MSTATI,13).EQ.0)  IACTVI=6
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.LT.0)             GO TO 27
      TEXT = 'File'
      IF (IACTVI.EQ.6)  TEXT='Data'
      WRITE (IQLOG,9024) LUNI,N,TEXT
 9024 FORMAT (1X/' FZIN.    LUN=',I4,' System EOF #',I4,
     F' seen as End-of-',A)
#endif
   27 LRTYP     = 0
      IQUEST(1) = IACTVI - 1
      GO TO 91

C----              Normal S/E-OF-RUN

   31 JRUNCR = LQ(LUHEAI)
      IF (JRUNCR.GT.0)             GO TO 34

C--                End of run

      IACTVI = MIN (4, MAX(IACTVI,2)+1)
      IF (JRUNCR.EQ.-1)  IACTVI = 4
      JRUNCR = IQ(KQSP+LQFI+29)
      IQUEST(1)  = IACTVI - 1
      IQUEST(11) = JRUNCR
      IQ(KQSP+LQFI+16-IACTVI) = IQ(KQSP+LQFI+16-IACTVI) + 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.LT.0)             GO TO 91
      IF (IACTVI.EQ.3)  THEN
          WRITE (IQLOG,9031) LUNI,JRUNCR
        ELSE
          WRITE (IQLOG,9032) LUNI
        ENDIF
 9031 FORMAT (1X/' FZIN.    LUN=',I4,' End of Run ',I8)
 9032 FORMAT (   ' FZIN.    LUN=',I4,' Zebra EoF')
#endif
      GO TO 91

C--                Start of run

   34 IQ(KQSP+LQFI+14) = IQ(KQSP+LQFI+14) + 1
      IQ(KQSP+LQFI+29) = JRUNCR
      IQUEST(1)  = 1
      IQUEST(11) = JRUNCR
      IACTVI = 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.0)  WRITE (IQLOG,9034) LUNI,JRUNCR
 9034 FORMAT (1X/' FZIN.    LUN=',I4,' Start of Run ',I6)
#endif
      GO TO 91

C-----------------------------------------------------
C-                 Other error conditions
C-----------------------------------------------------

   40 GO TO ( 41, 44, 64, 61, 71, 74, 77, 67), JRETCD

C--                Unexpected end-of-file

   41 IQ(KQSP+LQFI+30) = 1
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.-2)  WRITE (IQLOG,9041) LUNI
 9041 FORMAT (1X/' FZIDIA.  LUN=',I4,'  Unexpected EoF')
#endif
      GO TO 71

C------            Unexpected start/end of run

   44 IQ(KQSP+LQFI+30) = 2
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.-2)  WRITE (IQLOG,9044) LUNI
 9044 FORMAT (1X/' FZIDIA.  LUN=',I4,'  Unexpected Start/End-of-Run.')
#endif
      GO TO 71

C----              User error

   61 IQUEST(11)= -1
      ERRMSG = 'BAD CALLING'
      GO TO 81

C----              NOT ENOUGH SPACE

   64 IQUEST(11)= -2
      ERRMSG = 'NOT ENOUGH SPACE'
      GO TO 80

C----              EMERGENCY STOP

   67 IQUEST(11)= -3
      ERRMSG = 'EMERGENCY STOP'
      GO TO 80

C----              BAD DATA

   71 IQUEST(11)= -3
      ERRMSG = 'BAD DATA'
      GO TO 80

C----              BAD CONSTRUCTION

   74 IQUEST(11)= -4
      ERRMSG = 'BAD CONSTRUCTION'
      GO TO 80

C----              READ ERROR

   77 IQUEST(11)= -5
      ERRMSG = 'READ ERROR'

C------            Print error message

   80 IQ(KQSP+LQFI+18) = IQ(KQSP+LQFI+18) + 1
   81 IQUEST(1) = IQUEST(11)
      IQUEST(3) = IQ(KQSP+LQFI+22)
      IQUEST(12)= JERROR
      IQUEST(13)= LRTYP
      LRTYP = 0
      LAST  = IQ(KQSP+LQFI+26)
      IF (LAST.LT.0)       IQUEST(1)=MIN(-6,LAST-1)
#if defined(CERNLIB_QPRINT)
      IF (LOGLVI.GE.-2)
     +  WRITE (IQLOG,9081,ERR=82) LUNI,ERRMSG
     +,                    IQUEST(1),NRECAI,IQUEST(3)
     +,                    (IQUEST(J+10),J=1,NWERR+3)
   82 CONTINUE

 9081 FORMAT (1X/' FZIDIA.  LUN=',I4,1X,A/
     F 10X,'IQUEST(1/2/3)=',I4,2I8/
     F 10X,'IQUEST(11-) = ',3I4,2I9,
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))
     F 2O23)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))
     F 2Z17)
#endif
      IF (LAST.LT.-50)
     +         CALL ZFATAM ('FZIDIA - too many consecutive errors.')

   91 CONTINUE
#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
