*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZALFA)
      SUBROUTINE FZOALN

C-    Write 80 characters from the card buffer, left-adjust rest
C-    service routine to FZOASC

C-    Flush the buffer if IFLEND .GT. 0

#include "zebra/zbcd.inc"
#include "zebra/zkrakcc.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsysh.inc"
C--------------    End CDE                             --------------

      CHARACTER    LINEW*80
      EQUIVALENCE (LINEW,CQHOLK)

      EQUIVALENCE  (LUN ,IQUEST(90))
      EQUIVALENCE  (JPUT,IQUEST(91)), (IFLEND,IQUEST(92))


      NU = 80
      IF (IFLEND.NE.0)             GO TO 41

C--                Translate to characters

   12 DO 14  J=2,NU
   14 IQCETK(J) = IQLETT(IQCETK(J)+1)

      CALL UH1TOC (IQCETK,LINEW,80)

      WRITE (LUN,'(A)',IOSTAT=ISW) LINEW
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 41
          IQUEST(3) = ISW
          IQUEST(4) = LUN
          CALL ZFATAL
        ENDIF
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.NE.0)  WRITE (IQLOG,9817) LINEW
 9817 FORMAT (' DEVZE FZOALN, send line=',A)
#endif

C--                Copy excess numbers to start of line

      IQCETK(1) = IQBLAN

      JPUT = JPUT - 80
      IF (JPUT.LE.0)  THEN
          JPUT = 1
          RETURN
        ELSE
          DO 24  J=1,JPUT
   24     IQCETK(J+1) = IQCETK(J+80)
          JPUT = JPUT + 1
          IF (IFLEND.EQ.0)   RETURN
        ENDIF

C----              End-of-record, flush buffer

   41 IF (JPUT.GT.80)              GO TO 12
      IF (IFLEND.LT.0)       RETURN
      IQCETK(1) = IQLESS
      IF (JPUT.EQ.80)              GO TO 12
      CALL VBLANK (IQCETK(JPUT+1),80-JPUT)
      NU = JPUT
      GO TO 12
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
