*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVFPRI

C-    Print pending header lines for ZVERIF

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/zvfaut.inc"
#include "zebra/zvfwkc.inc"
C--------------    END CDE                             --------------
      DIMENSION    MBKDES(9)
      EQUIVALENCE (MBKDES(1),LBKLS)
      CHARACTER    CHWK*12

C--                Print main header line

      IF (IFLHD.EQ.0)  THEN
          IFLHD = 7
          N = LNBLNK (TEXTHD)
          WRITE (IQLOG,9011) TEXTHD(1:N)
          IF (IQVSTA.NE.0) WRITE (IQLOG,9012) IQVID
        ENDIF

 9011 FORMAT (1X/' ZVERIF.  called, ',A)
 9012 FORMAT (10X,'Verify ID =',I4,I10)

      IF (JVMODE.LT.0)             RETURN

C--                Print store header line

      IF (IFLST.EQ.0)  THEN
          IFLST = 7
          WRITE (IQLOG,9014) JQSTOR,(NQSNAM(J),J=1,2)
        ENDIF

 9014 FORMAT (1X/' ZVERIF.  Doing store',I3,', id= ',2A4)

      IF (JVMODE.EQ.0)             RETURN

C----              Link areas / banks

      LFAILS = JQSTOR
      LFAILD = JQDIVI
      IF (JVMODE.GE.2)             GO TO 31

C--                Doing link area

      IF (IFLLA.NE.0)              RETURN
      IFLLA = 7
      IF (LATMP.EQ.0)  THEN
          CHWK = ', permanent '
        ELSE
          CHWK = ', temporary '
        ENDIF
      WRITE (IQLOG,9024) JQSTOR,CHWK,NAMEPR,JDESLA
      RETURN

 9024 FORMAT (' ZVERIF.  Store',I3,A,'link-area ',2A4,' J=',I6)

C--                Doing banks

   31 IF (IFLDV.EQ.0)  THEN
          IFLDV = 7
          WRITE (IQLOG,9031) JQSTOR,JQDIVI,NAMEPR,LDIVLO,LDIVHI
        ENDIF

 9031 FORMAT (' ZVERIF.  Store',I3,', doing division',I3,1X,2A4,
     F' from/to =',2I10)

      IF (IFLBK.NE.0)              RETURN
      IFLBK = 7
      IF (IBKXX.EQ.0)  THEN
          CHWK(1:6) = ' live '
        ELSE
          CHWK(1:6) = ' dead '
        ENDIF
      WRITE (IQLOG,9034) CHWK(1:6),MBKDES
      RETURN

 9034 FORMAT (' ZVERIF.  Doing',A,'bank at',I9,' IDH/IDN= ',A4,I5,
     F' NIO=',I3/10X,'NL/NS/ND=',2I6,I8,' LA/LE=',2I10)
      END
*      ==================================================
#include "zebra/qcardl.inc"
