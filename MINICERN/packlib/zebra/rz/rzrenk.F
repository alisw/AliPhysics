*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:07  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZRENK(KEYOLD,KEYNEW)
*
************************************************************************
*
*        To rename a key in the CWD
* Input:
*   KEYOLD  Key array of dimension NWKEY containing the old key vector
*   KEYNEW  Key array of dimension NWKEY containing the new key vector
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 16.05.86
*  Last mod: 16.05.86
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
      DIMENSION KEYOLD(*),KEYNEW(*)
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
*           Check permission
*
      IFLAG=1
      CALL RZMODS('RZRENK',IFLAG)
      IF(IFLAG.NE.0)GO TO 99
*
      NKEYS=IQ(KQSP+LCDIR+KNKEYS)
      NWKEY=IQ(KQSP+LCDIR+KNWKEY)
      IF(NKEYS.LE.0)GO TO 90
*
      DO 10 K=1,NWKEY
         IKDES=(K-1)/10
         IKBIT1=3*K-30*IKDES-2
         IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
            KEY(K)=KEYOLD(K)
         ELSE
            CALL ZHTOI(KEYOLD(K),KEY(K),1)
         ENDIF
  10  CONTINUE
      DO 40 I=1,NKEYS
         DO 20 K=1,NWKEY
            LKC=IQ(KQSP+LCDIR+KLK)+(NWKEY+1)*(I-1)
            IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GO TO 40
  20     CONTINUE
         DO 30 K=1,NWKEY
            IKDES=(K-1)/10
            IKBIT1=3*K-30*IKDES-2
            IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
               IQ(KQSP+LCDIR+LKC+K)=KEYNEW(K)
            ELSE
               CALL ZHTOI(KEYNEW(K),IQ(KQSP+LCDIR+LKC+K),1)
            ENDIF
  30     CONTINUE
         GO TO 99
  40  CONTINUE
*
  90  IQUEST(1)=1
      IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,1000)
 1000 FORMAT(' RZRENK. Current directory contains no keys')
*
  99  RETURN
      END
