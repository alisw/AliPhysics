*
* $Id$
*
* $Log$
* Revision 1.5  1998/08/06 17:17:06  mclareni
* For multiple files, LUN has to be set in the common RZCLUN before it is used by RZFREE
*
* Revision 1.4  1997/07/04 15:26:01  couet
* - the error message printing previously added was wrong
*
* Revision 1.3  1997/07/03 09:21:15  couet
* - CHPATH is converted is uppercase before being compared to CHNAME
*
* Revision 1.2  1996/04/24 17:26:42  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZCLOS(CHPATH,CHOPT)
*
************************************************************************
*
*        To close all transactions with file CHPATH
*          Corresponding directories are dropped
*        A FORTRAN or CFCLOS is also issued for all associated files
* Input:
*   CHPATH  Character variable specifying the name of the top directory
*   CHOPT   Character variable specifying the options required
*       'A' Close all files currently open
*
* Called by <USER>
*
*  Author  : J. Shiers
*  Written : 11.11.91
*  Last mod: 11.11.91
*
************************************************************************
      CHARACTER*(*)     CHPATH,CHOPT
      CHARACTER*16      CHNAME,CHPATU
      DIMENSION         IHDIR(4)
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/zunit.inc"
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      LP = LENOCC(CHPATH)
      LC = LENOCC(CHOPT)
      IOPTA = 0
      IF(LC.GT.0) IOPTA = INDEX(CHOPT(1:LC),'A')
 
      IF(LQRS.EQ.0) RETURN
      LRZ=LQRS
  10  IF(LRZ.EQ.0)  RETURN
      LUN = IQ(KQSP+LRZ-5)
      IF(LUN.NE.0) THEN
         LOGLV = JBYT(IQ(KQSP+LRZ),15,3)-3
         CALL ZITOH(IQ(KQSP+LRZ+1),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHNAME,16)
         LN  = LENOCC(CHNAME)
*
*     Check top directory name unless IOPTA
*
         IF(IOPTA.EQ.0) THEN
            CHPATU = CHPATH(1:LP)
            CALL CLTOU(CHPATU)
            IF(CHPATU(1:LP).NE.CHNAME(1:LN)) GOTO 20
         ENDIF
         CALL RZEND(CHNAME(1:LN))
*
*    Close
*
         IF(LUN.GT.0) THEN
            IF(JBIT(IQ(KQSP+LRZ),5).EQ.0) THEN
               IF(LOGLV.GT.0) WRITE(IQLOG,*) 'RZCLOS. close unit ',LUN,
     +            ' (FORTRAN)'
               CLOSE(LUN)
            ELSE
#if defined(CERNLIB_QMIBM)
               IF(LOGLV.GT.0) WRITE(IQLOG,*) 'RZCLOS. C I/O not ',
     +            'supported on this system'
#endif
#if !defined(CERNLIB_QMIBM)
               IF(LOGLV.GT.0) WRITE(IQLOG,*) 'RZCLOS. close unit ',LUN,
     +            ' (C)'
               CALL CFCLOS(LUN-1000,0)
#endif
            ENDIF
         ENDIF
      ENDIF
  20  CONTINUE
      LRZ=LQ(KQSP+LRZ)
      GO TO 10
      END
