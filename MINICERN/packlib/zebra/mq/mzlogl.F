*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:41  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZLOGL (IXSTOR,LEVP)

C-    Set MZ logging level, user called

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------
      DIMENSION    LEVP(9)

#include "zebra/q_jbyt.inc"

      LEV = LEVP(1)
      IF (LEV.GE.100)              GO TO 21

#include "zebra/qstore.inc"

      LEV = MAX (-3,LEV)
      LEV = MIN  (4,LEV)
      LEP = MAX (LEV,NQLOGL)

      NQLOGM = LEV
      NQLOGL = LEV
      IQTABV(KQT+10) = LEV
#if defined(CERNLIB_QDEBPRI)
      IF (LEP.GE.0)  WRITE (IQLOG,9019) LEV,JQSTOR
 9019 FORMAT (1X/' MZLOGL.  Set Log Level',I2,' for store',I3)
#endif
      RETURN

   21 CONTINUE
#if defined(CERNLIB_QDEVZE)
      NQDEVZ = LEV - 100
      NQLOGM = MIN(NQDEVZ,2)
      NQLOGL = MIN(NQDEVZ,2)
#endif
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
