*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:28  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZICHH (IOMODE,MREC,IHOW)

C-    Check physical block header in MREC

C-     in :  IHOW = 1  fast block expected
C-                  0  steering block expected
C-                 -1  scanning for next steering block
C-                 -2  steering block expected when starting

C-    out :  IRETN in IQUEST(1)
C-           IRETN = 0  as expected
C-                   1  trouble
C-                   3  steering block ok, but with wrong length
C-                   4  end-of-file reached (only IBM just now)


#include "zebra/zunit.inc"
#include "zebra/quest.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzci.inc"
#include "fzhci.inc"
C--------------    End CDE                             --------------
      DIMENSION    MREC(8)
#include "fzstamp.inc"
#include "fzntolds.inc"
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBMD)
C     PARAMETER   (MASK   = Z FF000000)
      PARAMETER   (MASK   =  -16777216)
#endif

#include "zebra/q_jbyt.inc"


      IRETN = 0
      CALL UCOPY (MREC,MPHRI,8)
      JMODI = JBYT(MPHRI(5),30,3)
      NWRI  = JBYT(MPHRI(5),1,28)
      IF (MPHRI(1) .NE. MCCW1)     GO TO 71
      IF (MPHRI(2) .NE. MCCW2)     GO TO 71
      IF (MPHRI(3) .NE. MCCW3)     GO TO 71
      IF (MPHRI(4) .NE. MCCW4)     GO TO 71
      IF (JMODI.GE.5)              GO TO 71
      IF (JMODI.EQ.3)              GO TO 71
      IF (NWRI.GE.NTOLDS)          GO TO 71

      IF (JRECI.LT.0)              GO TO 71
      IF (NTLRI.LT.0)              GO TO 71
      IF (NTLRI.GT.NWRI)           GO TO 71
      IF (NFSTI.LT.0)              GO TO 71
      IF (NFSTI.GE.16384)          GO TO 71
      IF (IOMODE.NE.0)             GO TO 26
      IF (NFSTI*MAXREI.GT.NTOLDS)  GO TO 71
      IF (NWRI.NE.MAXREI)          GO TO 61

   26 IF (IHOW.GT.0)               GO TO 41
   27 IQUEST(1) = IRETN
      RETURN

C----              Unexpected steering block

   41 IRETN = 1
      GO TO 27

C----              Steering block of unexpected length

   61 IF (IHOW.LE.0)  IRETN = 3
      GO TO 27

C----              Not a steering block

   71 IF (IHOW.GT.0)               GO TO 27
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMIBMD)

C----     IBM :  detect end-of-file on direct-access file if :
C--                1) file-format D
C--                2) just after Zebra EoF
C--                3) record starts with FF000000, rest zero

      IF (IOMODE.NE.0)             GO TO 74
      IF (IFIFOI.NE.2)             GO TO 74
      IF (IACTVI.NE.4)             GO TO 74
      IF (MPHRI(1).NE.MASK)        GO TO 74
      DO 73  J=2,8
      IF (MREC(J).NE.0)            GO TO 74
   73 CONTINUE
      IRETN = 4
      GO TO 27

#endif
   74 IRETN = 1
#if defined(CERNLIB_QPRINT)
      IF (IHOW.EQ.-1)              GO TO 27
      IF (LOGLVI.LT.-1)            GO TO 27
      LUN = LUNI
      IF (IOMODE.NE.0)  LUN = LUNX
      IQUEST(1) = MCCW1
      IQUEST(2) = MCCW2
      IQUEST(3) = MCCW3
      IQUEST(4) = MCCW4
      CALL UCOPY (MPHRI,  IQUEST(5), 8)
      IQUEST(13) = NWRI
      CALL UCOPY (JRECI,  IQUEST(14),3)
      WRITE (IQLOG,9801) LUN,(IQUEST(J),J=1,16)
 9801 FORMAT (1X/' FZIPHx.  LUN=',I4,', Error 201, dump Block Header :'
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))&&(defined(CERNLIB_B60M))
     F/4X,'Must be:',4O23
     F/4X,'Oct 1-4:',4O23/8X,'5-8:',4O23/4X,'Dec 5-8:',4I23)
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))&&(!defined(CERNLIB_B60M))
     F/4X,'Must be:',4O17
     F/4X,'Oct 1-4:',4O17/8X,'5-8:',4O17/4X,'Dec 5-8:',4I17)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))
     F/4X,'Must be:',4Z17
     F/4X,'Hex 1-4:',4Z17/8X,'5-8:',4Z17/4X,'Dec 5-8:',4I17)
#endif
      GO TO 27
      END
