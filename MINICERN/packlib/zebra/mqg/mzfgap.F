*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:30  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZFGAP

C-    Find the 4 largest gaps,
C-    gaps 3 + 4 may be inside the moving memory region
C-                   defined by JQSTMV,JQDVM1+2,NQDVMV
C-    gaps 1 + 2 must be outside the moving memory region

C-    Entry 'J' in the table contains :
C-      IQGAP(1,J) = number of words available
C-      IQGAP(2,J) = abs. zebra adr of first word in the gap
C-      IQGAP(3,J) = division number
C-      IQGAP(4,J) = store number
C-      IQGAP(5,J)   user flag

C-    Original meaning of JQSTMV,JQDVM1,JQDVM2,NQDVMV (cf. MZTABS):
C-    unless JQSTMV<0, the divisions JQDVM1 to JQDVM2 inclusive
C-    of the store JQSTMV may have to be shifted by NQDVMV words,
C-    to the left if -ve, to the right if +ve.
C-    NQDVMV=0 means to the left by a yet unknown amount.
C-    It follows that the gaps after divisions JDIV are not available
C-    for the relocation table because they may be over-written
C-    by the memory move, with :
C-       NQDVMV  0, -ve :  JDIV = <JQDVM1-1,JQDVM2-1>
C-                  +ve :  JDIV = <JQDVM1,  JQDVM2  >

C-    MZFGAP is also used for I/O to find a gap where to put the
C-    relocation table on input. In this case there is just one
C-    gap unavailable, namely the one which will receive the data.
C-    This is the gap before a reverse division (JQMODE=1) or after
C-    a forward division (JQMODE=0). To block this gap one can give
C-    JQDVM1=JQDVM2=JDIV and NQDVMV=1-2*JQMODE, setting JQSTMV also.

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    NGAPV(7), JDIVV(7), JSTOV(7), JPV(7)
      EQUIVALENCE (NGAPV(1),IQUEST(11)), (JDIVV(1),IQUEST(21))
      EQUIVALENCE (JSTOV(1),IQUEST(31)), (JPV(1),  IQUEST(41))
      DIMENSION    NQGAPV(2)
      EQUIVALENCE (NQGAPV(1),NQGAPN)


      CALL VZERO (IQGAP,20)

   15 DO 17  J=1,6
   17 NGAPV(J) = 0

      IF (JQSTMV.LT.0)             GO TO 19
      KT     = NQOFFT(JQSTMV+1)
      JDVSH1 = JQDVM1
      JDVSH2 = JQDVM2
      IF (NQDVMV.GT.0)             GO TO 19
      IF (JDVSH1.EQ.IQTABV(KT+9))  JDVSH1=IQTABV(KT+8)+1
      JDVSH1 = JDVSH1 - 1
      JDVSH2 = JDVSH2 - 1
   19 MINGN  = 0
      MINGV  = 0
      JMINGN = 1
      JMINGV = 5

C------            Loop all stores

      JSTO = -1
   21 JSTO = JSTO + 1
      IF (JSTO.GT.NQSTOR)          GO TO 61
      IF (NQALLO(JSTO+1).LT.0)     GO TO 21

      KT   = NQOFFT(JSTO+1)
      JDVN = 1

C----              Loop all divisions

   31 JDIV = JDVN
      IF (JDIV.EQ.21)              GO TO 21
      JDVN = JDIV + 1
      IF (JDIV.EQ.IQTABV(KT+8))  JDVN=IQTABV(KT+9)

      NWGAP = LQSTA(KT+JDVN) - LQEND(KT+JDIV)
      IF (NWGAP.LT.164)            GO TO 31
      IF (JSTO.NE.JQSTMV)          GO TO 41
      IF (JDIV.LT.JDVSH1)          GO TO 41
      IF (JDIV.GT.JDVSH2)          GO TO 41

C--                Gap in moving region

      IF (NWGAP.LE.MINGV)          GO TO 31
      NGAPV(JMINGV) = NWGAP
      JDIVV(JMINGV) = JDIV
      JSTOV(JMINGV) = JSTO
      JMINGV = 5
      MINGV  = NGAPV(5)
      IF (MINGV.LE.NGAPV(6))       GO TO 31
      JMINGV = 6
      MINGV  = NGAPV(6)
      GO TO 31

C--                Gap in non-moving region

   41 IF (NWGAP.LE.MINGN)          GO TO 31
      NGAPV(JMINGN) = NWGAP
      JDIVV(JMINGN) = JDIV
      JSTOV(JMINGN) = JSTO
      JMINGN = 1
      MINGN  = NGAPV(1)
      DO 44  J=2,4
      IF (MINGN.LE.NGAPV(J))       GO TO 44
      JMINGN = J
      MINGN  = NGAPV(J)
   44 CONTINUE
      GO TO 31

C-------           Store the result

   61 DO 62  J=1,6
   62 JPV(J) = J

C--                Sort non-moving gaps

      JG = 1
   65 JF = JPV(JG)
      JN = JPV(JG+1)
      IF (NGAPV(JF).LT.NGAPV(JN))  GO TO 67
      IF (JG.EQ.3)                 GO TO 71
   66 JG = JG + 1
      GO TO 65

   67 JPV(JG)   = JN
      JPV(JG+1) = JF
      IF (JG.EQ.1)                 GO TO 66
      JG = JG - 1
      GO TO 65

C--                Sort moving gaps

   71 JG = 4
   75 JF = JPV(JG)
      JN = JPV(JG+1)
      IF (NGAPV(JF).LT.NGAPV(JN))  GO TO 77
      IF (JG.EQ.5)                 GO TO 81
   76 JG = JG + 1
      GO TO 75

   77 JPV(JG)   = JN
      JPV(JG+1) = JF
      IF (JG.EQ.3)                 GO TO 76
      JG = JG - 1
      GO TO 75

C--                Store 4 biggest gaps

   81 NQGAPN = 0
      NQGAP  = 0
      JSEL   = 1
      DO 87  JG=1,4
      JU   = JPV(JG)
      NWGAP= NGAPV(JU)
      IF (NWGAP.EQ.0)              GO TO 87
      JDIV = JDIVV(JU)
      JSTO = JSTOV(JU)
      KT   = NQOFFT(JSTO+1)
      KS   = NQOFFS(JSTO+1)
      IQGAP(1,JG) = NWGAP
      IQGAP(2,JG) = KS+ LQEND(KT+JDIV)
      IQGAP(3,JG) = JDIV
      IQGAP(4,JG) = JSTO
      IF (JU.GE.5)  JSEL=2
      NQGAPV(JSEL) = JG
   87 CONTINUE
      NQGAP = MAX (NQGAPN,NQGAP)
#if defined(CERNLIB_QDEVZE)
      IF (NQDEVZ.GE.11)
     +WRITE (IQLOG,9888) NQGAPV,(JG,(IQGAP(J,JG),J=1,5),JG=1,NQGAP)
 9888 FORMAT (1X/' DEVZE MZFGAP,  NQGAPN,NQGAP=',2I4
     F/16X,' JGAP    SIZE     LOC JDIV JSTO USED'
     F/(16X,I4,2I8,3I5))
#endif
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
