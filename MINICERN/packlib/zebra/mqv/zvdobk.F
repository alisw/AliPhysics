*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVDOBK

C-    Subsidiary to ZVERIF, check the links of one bank

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/zvfwkc.inc"
C--------------    End CDE                             --------------
      EQUIVALENCE (NL,NBKNN(2)), (NS,NBKNN(3))
      CHARACTER    CHWK*24, CHYN(4)*1

#include "zebra/q_jbit.inc"
      DATA  CHYN   / 'n', 'y', 'i', '.' /

#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  CALL ZVFPRI
#endif

C----              Check all links

      JLNK = -2
   22 JLNK = JLNK + 1
      IF (JLNK.GT.NL)              GO TO 61
      LINK = LQ(KQS+LBKLS-JLNK)
      IF (LINK.EQ.0)               GO TO 22
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  WRITE (IQLOG,9823) JLNK,LINK
 9823 FORMAT (16X,'Link number',I6,' is=',I10)
#endif
      IF (LINK.LT.LSTOLO)          GO TO 37
      IF (LINK.GE.LSTOHI)          GO TO 37

      IF (LINK.GE.LDIVLO .AND. LINK.LT.LDIVHI)  GO TO 24
      JDV = MZFDIV (-7,LINK)
      IF (JDV.EQ.0)                GO TO 31
      IF (JLNK.EQ.0)               GO TO 33
      IF (JBIT(MASKTO,JDV).EQ.0)   GO TO 41

   24 IF (JLNK.GT.NS)              GO TO 22

C--                Check structural link

      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.EQ.0)             GO TO 22

      NFATAL = NFATAL + 1
      CALL ZVFPRI
      IF (JLNK.NE.-1)  THEN
          WRITE (IQLOG,9027) CHWARN,JLNK,LINK
        ELSE
          WRITE (IQLOG,9028) CHWARN,LINK
        ENDIF
      GO TO 22

 9027 FORMAT (A,'Structural link',I6,' of value',I10,
     F' does not point to a bank')
 9028 FORMAT (A,'The "up" link of value',I10,
     F' does not point to a bank')

C--                Link points into a reserve area

   31 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9031) CHWARN,JLNK,LINK
 9031 FORMAT (A,'Link',I6,' of value',I10,
     F' does not point into a division')
      GO TO 22

C--                'next' link points to a different division

   33 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9033) CHWARN,LINK
 9033 FORMAT (A,'The "next" link of value',I10,
     F' points out of the division')
      GO TO 22

C--                Link points outside the bank region

   37 IF (NQOPTS(1).GE.2)          GO TO 22
      CALL ZVFPRI
      WRITE (IQLOG,9037) CHWARN,JLNK,LINK
 9037 FORMAT (A,'Link',I6,' of value',I10,
     F' does not point into the bank region')

      IF (NQOPTS(1).EQ.0)  THEN
          NFATAL = NFATAL + 1
        ELSE
          NWARN = NWARN + 1
        ENDIF
      GO TO 22

C--                Link points into non-reference division

   41 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9041) CHWARN,JLNK,LINK,JDV
 9041 FORMAT (A,'Link',I6,' of value',I10,' points to division',I3
     F/10X,'which is not declared in the cross-reference table,')

      CALL MZXRPR (-7)
      GO TO 22

C------            Check "up" and "origin" of a live bank

   61 IF (IBKXX.NE.0)              RETURN
      LUP  = LQ(KQS+LBKLS+1)
      LORG = LQ(KQS+LBKLS+2)
      IF (LORG.EQ.0)               GO TO 71
      IF (LORG.LT.IQTABV(KQT+13))     GO TO 64
      IF (LORG.GT.IQTABV(KQT+14))     GO TO 64
      IF (LQ(KQS+LORG).NE.LBKLS)      GO TO 67
      RETURN

C--                Origin link illegal

   64 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9064) CHWARN,LORG
      RETURN

 9064 FORMAT (A,'Origin link adr=',Z10,'x is out-of-bounds')

C--                LQ(LORG) does not point to the bank

   67 NFATAL = NFATAL + 1
      CALL ZVFPRI
      J = LQ(KQS+LORG)
      WRITE (IQLOG,9067) CHWARN,LORG,J,LBKLS
      RETURN

 9067 FORMAT (A,'Origin link adr=',Z10 ,'x'
     F/10X,'value LQ(LORG) =',Z10,'x is not LBK=',Z10,'x')

C--                LORG is zero

   71 IF (LUP.EQ.0)                RETURN
      NWARN = NWARN + 1
      CALL ZVFPRI
      WRITE (IQLOG,9071) CHWARN,LUP
      RETURN

 9071 FORMAT (A,'The "origin" link is zero, but the "up" link is',I10)
      END
*      ==================================================
#include "zebra/qcardl.inc"
