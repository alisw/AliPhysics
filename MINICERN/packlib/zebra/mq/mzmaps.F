*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:42  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZMAPS

C-    Print the map of divisions in the current store

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"


      IF (JQSTOR.LT.0)             RETURN
      IF (JQSTOR.GT.NQSTOR)        RETURN

      WRITE (IQLOG,9021) JQSTOR,NQSNAM(1),NQSNAM(2)
 9021 FORMAT (1X/'   Map of store',I3,1X,2A4/3X,24('-')/
     F1X/3X,'Division',40X,'Number of times'
     F/14X,'Kind',16X,'Max-size',16X,'Garb-coll.'
     F/12X,'Mode    Position      used   allowed',
     F3X,'Wiped  user  auto Pushd Redcd'/1X)

      JDIV = 1
   32 IF (JDIV.EQ.JQDVLL+1)  JDIV=JQDVSY

      MODE = JBIT (IQMODE(KQT+JDIV), 1)
      KIND = JBYT (IQKIND(KQT+JDIV),21,4)
      LSTA = LQSTA(KQT+JDIV)
      LEND = LQEND(KQT+JDIV)
      NINI = NQDINI(KQT+JDIV)
      NMAX = NQDMAX(KQT+JDIV)
      NAM1 = IQDN1(KQT+JDIV)
      NAM2 = IQDN2(KQT+JDIV)

      NSIZ = NQDSIZ(KQT+JDIV)
      NWIP = NQDWIP(KQT+JDIV)
      NGAU = NQDGAU(KQT+JDIV)
      NGAF = NQDGAF(KQT+JDIV)
      NPSH = NQDPSH(KQT+JDIV)
      NRED = NQDRED(KQT+JDIV)

      NWIP = MIN(NWIP,9999999)
      NGAU = MIN(NGAU,99999)
      NGAF = MIN(NGAF,99999)
      NPSH = MIN(NPSH,99999)
      NRED = MIN(NRED,99999)

      NSIZ = MAX (NSIZ,LEND-LSTA)
      NQDSIZ(KQT+JDIV) = NSIZ

      JPOS = LSTA
      IF (MODE.NE.0) JPOS=LEND

      WRITE (IQLOG,9032) JDIV,NAM1,NAM2,MODE,KIND
     +,                  JPOS,NSIZ,NMAX
     +,                  NWIP,NGAU,NGAF,NPSH,NRED
 9032 FORMAT (3X,I2,1X,2A4,2I2,3I10,I8,4I6)

      JDIV = JDIV + 1
      IF (JDIV.LE.20)              GO TO 32
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
