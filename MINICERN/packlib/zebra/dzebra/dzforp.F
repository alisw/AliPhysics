*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZFORP
      SAVE KFOTYP
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzioc.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bkfoparq.inc"

      CHARACTER   CHROUT*(*),CHSTAK*6,KFOTYP(0:11)*1
      PARAMETER (CHROUT = 'DZFORP')
      DATA KFOTYP /'U','B','I','F','D','H','*','S','*','N','*','L'/
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif





      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IPOS = 38

      DO 100 JFOCUR = 1,JFOEND,2


          IF (JFOCUR.EQ.JFOREP+1) THEN
              WRITE(CQLINE(IPOS:IPOS+1),'(''/ '')')
              IPOS = IPOS + 2
          ENDIF
          ITYPE  = MFO(JFOCUR)
          IF (ITYPE.EQ.IFOSEQ) THEN
              WRITE(CQLINE(IPOS:IPOS+1),'(''*'',A1)') KFOTYP(ITYPE)
              IPOS = IPOS + 3
                                                           GO TO 100
          ENDIF
          NWSEC  = MFO(JFOCUR+1)
          ITYPE  = MIN(ITYPE,8)
          IF (NWSEC.LT.0) THEN
              WRITE(CQLINE(IPOS:IPOS+1),'(''-'',A1)') KFOTYP(ITYPE)
              IPOS = IPOS + 3
          ELSEIF (NWSEC.EQ.0) THEN
              WRITE(CQLINE(IPOS:IPOS+1),'(''*'',A1)') KFOTYP(ITYPE)
              IPOS = IPOS + 3
          ELSE
              DO 10 I=1,100
                  IF(NWSEC.EQ.0)                           GO TO 20
                      IQUEST(I)=MOD(NWSEC,10)
                      NWSEC    = NWSEC/10
   10         CONTINUE
   20         DO 30 J=1,I-1
                  WRITE(CQLINE(IPOS:IPOS),'(I1)') IQUEST(I-J)
                  IPOS = IPOS + 1
   30         CONTINUE
              WRITE(CQLINE(IPOS:IPOS),'(A1)') KFOTYP(ITYPE)
              IPOS = IPOS + 2
          ENDIF

          IF (IPOS.GT.100) THEN
              CALL DZTEXT(0,CDUMMQ,1)
              IPOS = 23
              CQLINE = ' '
          ENDIF
  100 CONTINUE

      IF (IPOS.GT.23) CALL DZTEXT(0,CDUMMQ,1)

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
