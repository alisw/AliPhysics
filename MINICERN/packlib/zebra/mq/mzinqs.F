*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:11:36  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZINQS (IXSTOR)

C- Return the attributes of division IXSTOR :
C-     IQUEST(1)  opmode    operation mode from NQALLO
C-     IQUEST(2)  NQSTRU    number of structural permanent links
C-     IQUEST(3)  NQREF     number of permanent links
C-     IQUEST(4)  NQLINK    number of permanent + working space links
C-     IQUEST(5)  LQSTA(1)  first word after the working space
C-     IQUEST(6)  LQSTA(21) first word after the store
C-     IQUEST(7)  NQMINR    minimum size of the reserve area after div. 1
C-     IQUEST(8)  LQ2END    minimum value of LQEND(2)
C-     IQUEST(9)  JQDVLL    last low division
C-    IQUEST(10)  JQDVSY    system division, also : first high division
C-    IQUEST(11)  NQLOGL    log level for this store
C-    IQUEST(12)  NQSNAM    chars 1:4 of the printing name of the store
C-    IQUEST(13)            chars 5:8
C-    IQUEST(14)  LOWORG    lowest  possible origin-adr for the store
C-    IQUEST(15)  HIGORG    highest possible origin-adr for the store
C-    IQUEST(16)  INALC     initial allocation to last division created
C-    IQUEST(17)  ITOFF     table-offset J of the table,
C-                          if contained at the end of the store

#include "zebra/zunit.inc"
#include "zebra/zvfaut.inc"
#include "zebra/mqsys.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZIN, 4HQS   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZINQS /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZINQS  ')
#endif

#include "zebra/q_jbyt.inc"


#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

#include "zebra/qstore.inc"

      CALL UCOPY (IQTABV(KQT+1), IQUEST(2), 16)
      IQUEST(1) = NQALLO(JQSTOR+1)
      IQUEST(2) = IQUEST(4)
      IQUEST(3) = IQUEST(5)
      IQUEST(4) = IQUEST(6)

      IQUEST(5) = LQSTA(KQT+1)
      IQUEST(6) = LQSTA(KQT+21)
#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.2)  WRITE (IQLOG,9071) JQSTOR,(IQUEST(J),J=1,17)

 9071 FORMAT (/' MZINQS.  Inquiry for store',I3,' returns :'
     F/10X,'opmode NQSTRU  NQREF'
     F,'   NQLINK    LSTA1   LSTA21   NQMINR   LQ2END'
     F/9X,3I7,5I9
     F/10X,'JQDVLL JQDVSY NQLOGL'
     F,'  name1/2   lowORG  highORG  initalc  toffset'
     F/9X,3I7,1X,2A4,4I9)
#endif
#include "zebra/qtrace99.inc"
      RETURN

C-  (3) NQSTRU    number of structural permanent links
C-  (4) NQREF     number of permanent links
C-  (5) NQLINK    number of permanent + working space links
C-      LQSTA(1)  first word after the working space
C-      LQSTA(21) first word after the store
C-  (6) NQMINR    minimum size of the reserve area after div. 1
C-  (7) LQ2END    minimum value of LQEND(2)
C-  (8) JQDVLL    last low division
C-  (9) JQDVSY    system division, also : first high division
C- (10) NQLOGL    log level for this store
C- (11) NQSNAM    chars 1:4 of the printing name of the store
C- (12)           chars 5:8
C- (13) LOWORG    lowest  possible origin-adr for the store
C- (14) HIGORG    highest possible origin-adr for the store
C- (15) INALC     initial allocation to last division created
C- (16) ITOFF     table-offset J of the table, if contained in store
      END
*      ==================================================
#include "zebra/qcardl.inc"
