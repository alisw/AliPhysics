*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:49  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZTABX

C-    Find all involved divisions, ie. all divisions which
C-    point into the moving divisions

#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------


#include "zebra/q_or.inc"
#include "zebra/q_jbytet.inc"

C--                Merge kinds of active divisions

      MERGE = 0
      LMT   = LQMTA

   22 IF (LQ(LMT+1).LT.2)          GO TO 27
      JDIV = LQ(LMT)
      MERGE  = IOR (MERGE, IQKIND(KQT+JDIV))
   27 LMT  = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 22

C--                Pick divisions pointing

      LMT   = LQMTA
   32 IF   (LQ(LMT+1))       38, 33, 37
   33 JDIV = LQ(LMT)
      IF (JBYTET(IQRTO(KQT+JDIV),MERGE,1,26).EQ.0)  GO TO 38
      IF (JBYTET(IQRNO(KQT+JDIV),MERGE,1,26).EQ.0)  GO TO 38
      LQ(LMT+1) = 1

   37 LQMTLU = LMT + 8
   38 LMT    = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 32
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
