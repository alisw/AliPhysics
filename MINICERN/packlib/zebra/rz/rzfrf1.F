*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:51  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZFRF1(LUNFZ,IOPTH)
*
************************************************************************
*
*        Copy into current directory from  a sequential FZ file
* Input:
*   LUNFZ   Logical unit number of the FZ sequential access file
*   IOPTH   1 read only the highest cycle from LUNFZ
*           0 read all cycles
*
* Called by <RZFRFZ>
*
*  Author  : R.Brun DD/US/PD
*  Written : 14.05.86
*  Last mod: 08.12.92 JDS. Change chopt to char*3
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*1 CHOPTA
      CHARACTER*3 CHOPT
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

*
*            Read next key
*
  10  NH=KNSIZE
      CALL FZIN(LUNFZ,JQPDVS,LRZ0,-1,' ',NH,IHEAD)
      IF(IQUEST(1).NE.0)GO TO 99
      IF(NH.EQ.1.AND.IHEAD(1).EQ.77)GO TO 99
      IF(IOPTH.NE.0.AND.IHEAD(2).NE.1)GO TO 10
      LFROM=LQ(KQSP+LRZ0-1)
      IBIT4=JBIT(IHEAD(3),4)
      IFORM=JBYT(IHEAD(3),1,3)
      NWKEY=IQ(KQSP+LCDIR+KNWKEY)
      DO 20 I=1,NWKEY
         IKDES=(I-1)/10
         IKBIT1=3*I-30*IKDES-2
         IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).GE.3)THEN
            CALL ZITOH(KEY(I),KEY(I),1)
         ENDIF
  20  CONTINUE
*
*            Data structure
*
      IF(IBIT4.NE.0)THEN
         CHOPTA='A'
      ELSE
         CHOPTA=' '
      ENDIF
      IF(IFORM.EQ.0)THEN
         CHOPT='LW'//CHOPTA
         CALL RZOUT(JQPDVS,LFROM,KEY,ICY,CHOPT)
      ELSE
*
*            Vector
*
         CHOPT=CHOPTA
         IF(IFORM.EQ.1)CHOPT='B'//CHOPTA
         IF(IFORM.EQ.2)CHOPT='I'//CHOPTA
         IF(IFORM.EQ.5)CHOPT='H'//CHOPTA
         NDATA=IQ(KQSP+LFROM-1)
         LTEMP=LQ(KQSP+LFROM)
         IF(LTEMP.NE.0)THEN
            NTOT=0
            NEW=1
  55        IF(NTOT.LT.NDATA)THEN
               CALL RZLINC(IQ(KQSP+LFROM+1),NTOT,IQ(KQSP+LTEMP+1),NEW)
               NEW=NEW+1
               GO TO 55
            ENDIF
         ENDIF
         CALL RZVOUT(IQ(KQSP+LFROM+1),NDATA,KEY,ICY,CHOPT)
      ENDIF
      IF(IQUEST(1).NE.0)THEN
        CALL MZDROP(JQPDVS,LFROM,'L')
        GO TO 99
      ENDIF
      LC=IQ(KQSP+LCDIR+KLC)
      IQ(KQSP+LCDIR+LC+KFLCYC)=IHEAD(3)
      CALL MZDROP(JQPDVS,LFROM,'L')
      GO TO 10
*
  99  RETURN
      END
