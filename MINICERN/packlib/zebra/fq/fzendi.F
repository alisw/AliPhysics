*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZENDI (LUNPAR,CHOPT)

C-    TERMINATION OF ZEBRA INPUT FILE

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzstat.inc"
#include "zebra/fzcf.inc"
#include "zebra/fzci.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNPAR(9)
      CHARACTER    CHOPT*(*)
      DIMENSION    IOPNUM(4)
      EQUIVALENCE (IOPTT,IOPTVF(1)), (IOPTR,IOPTVF(6))
     +,           (IOPTN,IOPTVF(2)), (IOPTU,IOPTVF(7))
     +,           (IOPTC,IOPTVF(3)), (IOPTX,IOPTVF(8))
     +,           (IOPTI,IOPTVF(4)), (IOPTK,IOPTVF(9))
     +,           (IOPTO,IOPTVF(5)), (IOPTQ,IOPTVF(10))
     +,           (IOPNUM(1),IOPTVF(11))
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZEN, 4HDI   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZENDI /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZENDI  ')
#endif

#include "zebra/q_jbit.inc"

#include "zebra/qtrace.inc"

      LUNORG  = LUNPAR(1)
      LUN     = LUNORG
      CALL UOPTC (CHOPT,'TNCIORUXKQ0123',IOPTVF)
      INFLUN = 0

      NEOF = -1
      IF (IOPNUM(4).NE.0)  NEOF=3
      IF (IOPNUM(3).NE.0)  NEOF=2
      IF (IOPNUM(2).NE.0)  NEOF=1
      IF (IOPNUM(1).NE.0)  NEOF=0
      IOPTK  = 1

C--                Terminate option selected

      IF (IOPTT.NE.0)  THEN
          IOPTN  = 0
          IOPTC  = 0
          IOPTI  = 0
          IOPTO  = 0
          IOPTK  = 0

C--                New option selected

        ELSEIF (IOPTN.NE.0)  THEN
          IOPTC  = 0
          IOPTI  = 0
          IOPTO  = 0
        ELSE
          IOPTR  = 0
          IOPTU  = 0
          IOPTX  = 0
        ENDIF

C--                Continue option selected

      IF (IOPTC.NE.0)  THEN
          IOPTI  = 0
          IOPTO  = 0
        ENDIF

C--                Input option selected

      IF (IOPTI.NE.0)  THEN
          IOPTO  = 0
          IOPTR  = 1
        ENDIF

C----              LOOP OVER ALL INPUT FILES IF LUNPAR=0

      IF (LUN.NE.0)                GO TO 24
      LQFF = LQFS
      GO TO 22

   21 LQFF = LQ(KQSP+LQFF)
   22 IF (LQFF.EQ.0)               GO TO 999
      LUN = IQ(KQSP+LQFF-5)

      LUNI = 0
      IQUEST(1) = -7
      CALL FZLOC (LUN,-1)
      IF (LUNI.EQ.LUN)             GO TO 31
      GO TO 21

C--                FZENDI called for one particular input file

   24 LUNI = 0
      IQUEST(1) = 0
      CALL FZLOC (LUN,-1)
      IF (LUNI.EQ.LUN)             GO TO 31
#if defined(CERNLIB_QPRINT)
      IF (IOPTQ.EQ.0)  THEN
          WRITE (IQLOG,9024) LUN
          IF (LUNF.EQ.LUN)  WRITE (IQLOG,9025) IACTVF
        ENDIF
 9024 FORMAT (1X/' FZENDI.  LUN=',I3,' is not a valid input file')
 9025 FORMAT (10X,'last activity=',I2)
#endif
      GO TO 999

C-------------------------------------------------
C-                 Do input file
C-------------------------------------------------

   31 IF (IOPTQ.NE.0)  LOGLVF=-2
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.0)
     +WRITE (IQLOG,9031) LUN,IACTVF,CHOPT
 9031 FORMAT (1X/' FZENDI.  For input file at LUN=',I3,
     F', Last activity=',I2,', OPT= ',A)
#endif
      IF (IACTVF.GE.8)             GO TO 61

C----              Print file usage

#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.LT.0)             GO TO 39
      IF (IQ(KQSP+LQFF+28).EQ.IQ(KQSP+LQFF+15))  GO TO 39
      N = 21
      IF (IFIFOF.NE.0)  N=24
      WRITE (IQLOG,9036)  (IQ(KQSP+LQFF+J),J=11,N)
 9036 FORMAT (10X,'Number of objects read : '
     F/10X,I9,' System EOF'
     F/10X,I9,' Zebra  EOF'
     F/10X,I9,' End-of-Run'
     F/10X,I9,' Start-of-Run'
     F/10X,I9,' Pilot records'
     F/10X,I9,' Non-empty d/s selected'
     F/10X,I9,' Empty d/s selected'
     F/10X,I9,' Read or Data errors'
     F/7X,I12,' Mega-words +'
     F/7X,I12,' words'
     F/10X,I9,' Good logical records',:/
     F 10X,I9,' Good physical records'
     F/10X,I9,' Steering blocks'
     F/10X,I9,' Words with conversion problems')

      IQ(KQSP+LQFF+28) = IQ(KQSP+LQFF+15)
#endif
   39 INFLUN = LUN
      INFSTA = IQ(KQSP+LQFF)
      CALL UCOPY (IQ(KQSP+LQFF+1), INFOFZ, 40)

C----              New file to be connected by the user

      IF (IOPTN.EQ.0)              GO TO 44
      IQ(KQSP+LQFF+2) = 8
      CALL VZERO (IQ(KQSP+LQFF+30),7)
      IF (IOPTR.NE.0)              GO TO 57
      GO TO 71

C--                Continue to read next file on same

   44 IF (IOPTC.EQ.0)              GO TO 51
      IF (IACTVF.NE.5)             GO TO 48
      IF (IACMOF.NE.1)             GO TO 79
      IF (JBIT(MSTATF,13).EQ.0)    GO TO 79

#include "fzendicx.inc"
* Ignoring t=pass
      GO TO 79

   48 CONTINUE
#if defined(CERNLIB_QPRINT)
      IF (LOGLVF.GE.-2)  WRITE (IQLOG,9048) LUNF
 9048 FORMAT (' FZENDI.  LUN=',I3,' not positioned on system EoF;',
     F' C-option is redundant.')
#endif
      GO TO 79

C----              REWIND

   51 IF (IOPTR.EQ.0)              GO TO 61

      IQ(KQSP+LQFF+2) = 8
      CALL VZERO (IQ(KQSP+LQFF+11),23)

   57 CALL FZMACH (0)

C----              Switch to output

   61 IF (IOPTO.EQ.0)              GO TO 71
      IQ(KQSP+LQFF+2) = 10
      IF (NEOF.GE.0)  THEN
          CALL SBYT (NEOF, IQ(KQSP+LQFF), 13,2)
#if defined(CERNLIB_QPRINT)
          IF (LOGLVF.GE.0)  WRITE (IQLOG,9062) LUN,NEOF
 9062 FORMAT (10X,'LUN=',I4,' Set NEOF =',I2)
#endif
        ENDIF

C----              Close the file, Drop the control bank

   71 IF (IOPTX.NE.0)  THEN
          CALL FZMACH (2)
          IQ(KQSP+LQFF+1) = IADOPF
        ENDIF

      IF (IOPTK.EQ.0)  CALL MZDROP (0,LQFF,'.')

   79 LUNI = 0
      IF (LUNORG.EQ.0)             GO TO 21

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
