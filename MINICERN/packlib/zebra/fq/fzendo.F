*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:20  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZENDO (LUNPAR,CHOPT)

C-    TERMINATION OF ZEBRA OUTPUT FILE

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzstat.inc"
#include "zebra/fzcf.inc"
#include "zebra/fzcx.inc"
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
      DATA  NAMESR / 4HFZEN, 4HDO   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZENDO /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZENDO  ')
#endif


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

      LEVOUT = 0
      IOPTK  = 1

C--                Terminate option selected

      IF (IOPTT.NE.0)  THEN
          LEVOUT = 16
          IOPTN  = 0
          IOPTC  = 0
          IOPTI  = 0
          IOPTO  = 0
          IOPTK  = 0

C--                New option selected

        ELSEIF (IOPTN.NE.0)  THEN
          LEVOUT = 16
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
          LEVOUT = 15
          IOPTI  = 0
          IOPTO  = 0
        ENDIF

C--                Input option selected

      IF (IOPTI.NE.0)  THEN
          LEVOUT = 16
          IOPTR  = 1
          IF (IOPTO.NE.0)  IOPTO=-1
        ENDIF

C--                Output option selected

      IF (IOPTO.GT.0)  THEN
          LEVOUT = 16
          IOPTR  = 1
        ENDIF

C----              LOOP OVER ALL OUTPUT FILES IF LUNPAR=0

      IF (LUN.NE.0)                GO TO 24
      LQFF = LQFS
      GO TO 22

   21 LQFF = LQ(KQSP+LQFF)
   22 IF (LQFF.EQ.0)               GO TO 999
      LUN = IQ(KQSP+LQFF-5)

      LUNX = 0
      IQUEST(1) = -7
      CALL FZLOC (LUN,-2)
      IF (LUNX.EQ.LUN)             GO TO 31
      GO TO 21

C--                FZENDO called for one particular output file

   24 LUNX = 0
      IQUEST(1) = 0
      CALL FZLOC (LUN,-2)
      IF (LUNX.EQ.LUN)             GO TO 31
#if defined(CERNLIB_QPRINT)
      IF (IOPTQ.EQ.0)   THEN
          WRITE (IQLOG,9024) LUN
          IF (LUNF.EQ.LUN)  WRITE (IQLOG,9025) IACTVF
        ENDIF
 9024 FORMAT (1X/' FZENDO.  LUN=',I3,' is not a valid output file')
 9025 FORMAT (10X,'last activity=',I2)
#endif
      GO TO 999

C-------------------------------------------------
C-                 Do output file
C-------------------------------------------------

   31 IF (IOPTQ.NE.0)  LOGLVX=-2
      LOGLVF = LOGLVX

C-    LEVOUT = 13   flush buffer
C-             14   write E-o-Run
C-             15   write E-o-File
C-             16   write E-o-Data

#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)
     +WRITE (IQLOG,9031) LUN,IACTVF,CHOPT
 9031 FORMAT (1X/' FZENDO.  For output file at LUN=',I3,
     F', Last activity=',I2,', OPT= ',A)
#endif

      IF (IACTVF.LT.11)            GO TO 71
      IF (IACTVF.EQ.18)            GO TO 71
      IF (MEDIUX.LT.6)  THEN
          IF (LEVOUT.GT.IACTVF)
     +            CALL FZOUT (LUN,-7,0,LEVOUT,'FZEND',0,0,0)
        ENDIF

C----              Print file usage

#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.LT.0)             GO TO 39
      IF (IQ(KQSP+LQFF+28).EQ.IQ(KQSP+LQFF+15))  GO TO 39
      N = 21
      IF (IFIFOF.NE.0)  N=24
      WRITE (IQLOG,9036)  (IQ(KQSP+LQFF+J),J=11,N)
 9036 FORMAT (10X,'Number of objects written : '
     F/10X,I9,' System EOF'
     F/10X,I9,' Zebra  EOF'
     F/10X,I9,' End-of-Run'
     F/10X,I9,' Start-of-Run'
     F/10X,I9,' Pilot records'
     F/10X,I9,' Non-empty d/s'
     F/10X,I9,' Empty d/s'
     F/10X,I9,' Number of errors'
     F/7X,I12,' Mega-words +'
     F/7X,I12,' words'
     F/10X,I9,' Logical records',:/
     F 10X,I9,' Physical records'
     F/10X,I9,' Steering blocks'
     F/10X,I9,' Words with conversion problems')

      IQ(KQSP+LQFF+28) = IQ(KQSP+LQFF+15)
#endif
   39 INFLUN = LUN
      INFSTA = IQ(KQSP+LQFF)
      CALL UCOPY (IQ(KQSP+LQFF+1), INFOFZ, 40)

C----              New file to be connected by the user

      IF (IOPTN.EQ.0)              GO TO 44
      IQ(KQSP+LQFF+2) = 10
      CALL VZERO (IQ(KQSP+LQFF+30),7)
      IF (IOPTR.NE.0)              GO TO 57
      GO TO 71

C----              Switch to Input

   44 IF (IOPTI.EQ.0)              GO TO 51
      IF (NEOF.GE.0)  THEN
          CALL SBYT (NEOF, IQ(KQSP+LQFF), 13,2)
#if defined(CERNLIB_QPRINT)
          IF (LOGLVX.GE.0)  WRITE (IQLOG,9044) LUN,NEOF
 9044 FORMAT (10X,'LUN=',I4,' Set NEOF =',I2)
#endif
        ENDIF

C----              Rewind

   51 IF (IOPTR.EQ.0)              GO TO 71

      CALL VZERO (IQ(KQSP+LQFF+11),23)

      IF (IOPTO.GT.0)  THEN
          IQ(KQSP+LQFF+2) = 10
        ELSE
          IQ(KQSP+LQFF+2) = 18
          IF (IOPTO.EQ.0)  CALL SBIT0 (IQ(KQSP+LQFF),12)
        ENDIF

   57 CALL FZMACH (IOPTU)

C----              Close the file, Drop the control bank

   71 IF (IOPTX.NE.0)  THEN
          CALL FZMACH (2)
          IQ(KQSP+LQFF+1) = IADOPF
        ENDIF

      IF (IOPTK.EQ.0)  CALL MZDROP (0,LQFF,'.')

      LUNX = 0
      IF (LUNORG.EQ.0)             GO TO 21

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
