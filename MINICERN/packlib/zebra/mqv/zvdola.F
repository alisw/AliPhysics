*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:13:26  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:13  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZVDOLA (MDESV)

C-    Subsidiary to ZVERIF, check all link areas

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/zvfwkc.inc"
C--------------    End CDE                             --------------
      DIMENSION    MDESV(99)

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


      JDESMX = MDESV(1) - 4
      JDESLA = -4
      IF (MDESV(2).GE.MDESV(3))   JDESLA =1

C------            Next link area

   17 JDESLA = JDESLA + 5
      IF (JDESLA.GE.JDESMX)        RETURN
      IF (JDESLA.EQ.6)             GO TO 17
      LOCAR  = MDESV(JDESLA+1)
      LOCARE = MDESV(JDESLA+2)
      MODAR  = MDESV(JDESLA+3)
      LATMP  = JBIT(MODAR,31)
      LOCARR = LOCAR + JBYT (MODAR,1,15)
      NAMEPR(1) = MDESV(JDESLA+4)
      NAMEPR(2) = MDESV(JDESLA+5)

      IFLLA = 0
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  THEN
          CALL ZVFPRI
          WRITE (IQLOG,9818) JDESLA,LOCAR,LOCARR,LOCARE,LATMP
        ENDIF
 9818 FORMAT (10X,'Area at',I5,' LA,LR,LE,temp=',3I10,I2)
#endif
      LIX = LOCAR
      IF (LATMP.EQ.0)              GO TO 21

C--       temporary link area

      LIX = LIX + 2
      JDMUST = LQ(KQS+LOCAR+1)
      IF (JDMUST.EQ.JDESLA)        GO TO 18
      NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9018) CHWARN,JDMUST,JDESLA
 9018 FORMAT (A,'Word 2 destroyed, is= ',Z8,'x, should be= ',Z8)

   18 IF (LQ(KQS+LOCAR).EQ.0)         GO TO 17

C----              Check all links

   21 LIX  = LIX - 1
      JLNK = 0

C--                Next link

   22 LIX = LIX + 1
      IF (LIX.EQ.LOCARE)           GO TO 17
      JLNK = JLNK + 1
      LINK = LQ(KQS+LIX)
      IF (LINK.EQ.0)               GO TO 22
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)  WRITE (IQLOG,9823) LINK,LIX
 9823 FORMAT (16X,'Link =',I9,' from LIX =',I9)
#endif
      IF (LINK.LT.LSTOLO)          GO TO 37
      IF (LINK.GE.LSTOHI)          GO TO 37

      JDV = MZFDIV (-7,LINK)
      IF (JDV.EQ.0)                GO TO 31
      IF (LIX.GE.LOCARR)           GO TO 22

C--                Check structural link

      CALL MZCHLS (-7,LINK)
      IF (IQFOUL.EQ.0)             GO TO 22

      NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9027) CHWARN,JLNK,LINK
 9027 FORMAT (A,'Structural link number',I6,' of value',I10,
     F' does not point to a bank')
      GO TO 22

C--                Link points into a reserve area

   31 NFATAL = NFATAL + 1
      CALL ZVFPRI
      WRITE (IQLOG,9031) CHWARN,JLNK,LINK
 9031 FORMAT (A,'Link number',I6,' of value',I10,
     F' does not point into a division')
      GO TO 22

C--                Link points outside the bank region

   37 IF (NQOPTS(1).GE.2)          GO TO 22
      CALL ZVFPRI
      WRITE (IQLOG,9037) CHWARN,JLNK,LINK
 9037 FORMAT (A,'Link number',I6,' of value',I10,
     F' does not point into the bank region')

      IF (NQOPTS(1).EQ.0)  THEN
          NFATAL = NFATAL + 1
        ELSE
          NWARN = NWARN + 1
        ENDIF
      GO TO 22
      END
*      ==================================================
#include "zebra/qcardl.inc"
