*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZIALN

C-    Read 80 characters to the card buffer;
C-    service routine to FZIASC.

#include "zebra/zbcd.inc"
#include "zebra/zkrakcc.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      CHARACTER    LINER*80
      EQUIVALENCE (LINER,CQHOLK)

      EQUIVALENCE   (LUN,IQUEST(90)), (JSKIP,IQUEST(91))
      EQUIVALENCE (MRSTA,IQUEST(92)), (MREND,IQUEST(93))
      EQUIVALENCE  (JTKC,IQUEST(96)),  (JTKL,IQUEST(97))
     +,            (JTKE,IQUEST(98))


C--                Save trailing card buffer not yet analysed

      JTKE  = JTKE - 1
      NSV   = JTKE - JTKC
      IWDSV = IQKRAK(JTKE)

      IF (NSV.GT.0)  THEN
#if defined(CERNLIB_FQABLANK)
          NSH = 160 - JTKE
          CALL UCOPY2 (IQKRAK(JTKC),IQKRAK(JTKC+NSH),NSV)
          JTKC = JTKC + NSH
          JTKE = 160
#endif
          DO 18  J=JTKC,JTKE-1
   18     IQKRAK(J-79) = IQKRAK(J)
        ENDIF

C--                Read next line

   21 CONTINUE
 8000 FORMAT (A)
#if defined(CERNLIB_QERREND)
      READ (LUN,8000,END=81,ERR=82) LINER
#endif
#if !defined(CERNLIB_QERREND)
      READ (LUN,8000,IOSTAT=JSTAT) LINER
      IF (JSTAT.NE.0)              GO TO 81
#endif
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.3)  WRITE (IQLOG,9821) LINER
 9821 FORMAT (' FZIALN-  Read line=',A)
#endif

      CALL IZCHAV (80)

      JSTA = 2
      NOFF = 6

C--                Check special conditions

      IF (IQCETK(1).NE.45)         GO TO 41
      IF (JSKIP.NE.0)              GO TO 51

C--                Ready start/end pointers

   26 IF (NSV.GE.0)  THEN
          JTKC = 81 - NSV
          IQCETK(1) = IWDSV
        ELSE
          JTKC = 80 + JSTA
        ENDIF
      JTKE = 161

C--                Purge blanks from line, if FQABLANK

#if defined(CERNLIB_FQABLANK)
      JTKE = 80 + JSTA
      DO 34  J=JSTA,80
      IF (IQCETK(J).EQ.45)         GO TO 34
      IQKRAK(JTKE) = IQCETK(J)
      JTKE = JTKE + 1
   34 CONTINUE
#endif
      JTKL = JTKE - NOFF
      RETURN

C----              Control-line read

C--       End of physical record

   41 IF (IQCETK(2).EQ.59)         GO TO 45
      IF (IQCETK(1).NE.60)         GO TO 91
      IF (JSKIP.NE.0)              GO TO 21
      MREND = 1
      NOFF  = 0
      GO TO 26

C--       Start of physical record

   45 NSV   = -1
      MRSTA = 1
      IF (IQCETK(3).NE.27)  MRSTA=2

      JSTA = 4
      IF (IQCETK(1).EQ.59)         GO TO 26
      IF (IQCETK(1).NE.60)         GO TO 91
      MREND = 1
      NOFF  = 0
      GO TO 26

C----              Fast skip of trailing record

   51 CONTINUE
#if defined(CERNLIB_QERREND)
      READ (LUN,8000,END=81,ERR=82) LINER(1:1)
#endif
#if !defined(CERNLIB_QERREND)
      READ (LUN,8000,IOSTAT=JSTAT) LINER(1:1)
      IF (JSTAT.NE.0)              GO TO 81
#endif
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVI.GE.4)  WRITE (IQLOG,9852)
 9852 FORMAT (' FZIALN-  Skip line')
#endif

      IF (LINER(1:1).EQ.' ')       GO TO 51
      GO TO 21

C------            End-of-File, Read Error

#if defined(CERNLIB_QERREND)
   81 IQUEST(1) = -2
   82 IQUEST(1) = IQUEST(1) + 1
#endif
#if !defined(CERNLIB_QERREND)
   81 IQUEST(1) = JSTAT
#endif
      RETURN

C----              Context error

   91 IQUEST(1) = 7799
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
