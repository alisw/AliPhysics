*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZIOTC (IXVECT,LBANK,NCHTR,CHTR)

C-    Reconstruct CHTR from the I/O decriptor :
C-        CALL MZIOTC (IXSTOR,L, NCHTR*, CHTR*) for bank at L
C-        CALL MZIOTC (IODESC,0, NCHTR*, CHTR*) for descriptor vector
C-    Return :
C-      NCHTR = useful length of CHTR, -ve if trouble
C-       CHTR = the text string

#include "zebra/mqsys.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzioc.inc"
C-------------------------  CDE  -----------------------------------
      COMMON /SLATE/ NDSLAT,DUMMY(39)
      DIMENSION    IXVECT(9)
      CHARACTER    CHTR*(*), CHINV*(8)
      DATA CHINV   / 'UBIFDHXS' /


C----              Crack the I/O desciptor

      IF (LBANK.EQ.0)  THEN
          CALL MZIOCR (IXVECT)
        ELSE
          CALL MZCHLS (IXVECT,LBANK)
          NT = 0
          IF (IQFOUL.NE.0)         GO TO 91
          CALL MZIOCR (LQ(KQS+IQLN))
        ENDIF

C----              Loop to convert cracked to printable
C--                (mark wanted blanks with !)

      NCHMX  = LEN(CHTR) - 2
      CHTR   = ' '
      JFOCUR = 0

      JT = 1
   21 IF (JFOCUR.EQ.JFOREP)  THEN
        CHTR(JT:JT+1) = '/!'
        JT = JT + 2
      ENDIF

      JF = MFO(JFOCUR+1)
      JC = MFO(JFOCUR+2)

      IF (JC.GT.0)  THEN
          NEED = 8
        ELSE
          NEED = 3
        ENDIF
      IF (JT-1+NEED.GT.NCHMX)      GO TO 61

C--                Store next sector code

   41 IF    (JC)             44, 45, 46
   44 CHTR(JT:JT) = '-'
      JT = JT + 1
      GO TO 48

   45 CHTR(JT:JT) = '*'
      JT = JT + 1
      GO TO 48

   46 CALL CSETDI (JC,CHTR,JT,JT+5)
      JT = JT + 6

   48 CHTR(JT:JT)= CHINV(JF+1:JF+1)
      CHTR(JT+1:JT+1) = '!'
      JT = JT + 2

      JFOCUR = JFOCUR + 2
      IF (JFOCUR.LT.JFOEND)        GO TO 21
      NEED = -1

C--                Squeeze unwanted blanks

   61 CALL CLEFT (CHTR,1,JT-1)
      NT = NDSLAT

      IF (NEED.GE.0)  THEN
          JT = NT + 1
          IF (NT+NEED.LE.NCHMX)    GO TO 41
        ENDIF

C--                Change to blank

      DO 63 JX=1,NT
      IF (CHTR(JX:JX).EQ.'!')  CHTR(JX:JX) = ' '
   63 CONTINUE

      IF (NEED.GE.0)  THEN
          CALL CLEFT (CHTR,1,NT)
          NT = NDSLAT
          JT = NT + 1
          IF (NT+NEED.LE.NCHMX)    GO TO 41
          NT = 0
        ENDIF

   91 NCHTR = NT
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
