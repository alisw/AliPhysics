*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:50  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZFREE(CHLOCK)
*
************************************************************************
*
*      To free a directory previously locked by RZLOCK
* Input:
*   CHLOCK  Character variable identifying the owner of the lock.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.05.86
*  Last mod: 09.01.91
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
      CHARACTER*(*) CHLOCK
      DIMENSION IHL(2)
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      IF(LTOP.EQ.0)GO TO 99
*
      LOGLV=JBYT(IQ(KQSP+LTOP),15,3)-3
      NCH=LEN(CHLOCK)
      IF(NCH.GT.8)NCH=8
      CALL UCTOH(CHLOCK,IHL,4,NCH)
      IF(NCH.LT.5)CALL VBLANK(IHL(2),1)
      CALL ZHTOI(IHL,IHL,2)
*
*           Check write permission
*
***      IF(JBIT(IQ(KQSP+LCDIR),1).NE.0)THEN
***         IQUEST(1)=4
***         IF(LOGLV.GE.-2) WRITE(IQLOG,9010)
*** 9010    FORMAT(' RZFREE. No authorisation to write in that directory')
***         GO TO 99
***      ENDIF
*
*           Lock first record
*
      LRIN  = LQ(KQSP+LTOP-7)
      LPURG = LQ(KQSP+LTOP-5)
      LROUT = LQ(KQSP+LTOP-6)
      IF(LRIN.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LRIN,LTOP,-7,'RZIN',0,0,LREC+1,2,-1)
         IQ(KQSP+LRIN-5)=IQ(KQSP+LTOP-5)
      ENDIF
      NWL =50
      NTRY=0
  10  CALL RZIODO(LUN,NWL,1,IQ(KQSP+LRIN+1),1)
      IF(IQUEST(1).NE.0)GO TO 99
      IF(IQ(KQSP+LRIN+2).GT.NWL)THEN
         NWL=IQ(KQSP+LRIN+2)
         GO TO 10
      ENDIF
      NWL=IQ(KQSP+LRIN+2)
      IQ(KQSP+LTOP+KIRIN)=0
      IF(IQ(KQSP+LRIN+3).NE.0)THEN
         NWL=50
         NTRY=NTRY+1
#if defined(CERNLIB_QMVAX)
         IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
         CALL LIB$WAIT(0.1)
#endif
         IF(NTRY.LT.100)GO TO 10
         IF(LOGLV.GE.-2) WRITE(IQLOG,1000)
 1000    FORMAT(' RZFREE. Cannot lock that directory')
         IQUEST(1)=1
         GO TO 99
      ENDIF
      IQ(KQSP+LRIN+3)=1
      CALL RZIODO(LUN,MAX(NWL,50),1,IQ(KQSP+LRIN+1),2)
*
*           Read fresh Top directory
*
      IF(LTOP.NE.LCDIR)CALL RZRTOP
*
*           Save directories
*
      CALL RZSAVE
*
      IQ(KQSP+LRIN+3)=0
      NLOCK=IQ(KQSP+LRIN+1)
      LDC=IQ(KQSP+LCDIR+KLD)
      IRD=IQ(KQSP+LCDIR+LDC+1)
*
*           Search lock-id
*
      IF(NLOCK.GT.0)THEN
         LL=4
  20     NLESS=IQ(KQSP+LRIN+LL)
         IF(NLESS.NE.0)THEN
            IF(IQ(KQSP+LRIN+LL+1).EQ.IHL(1).AND.
     +         IQ(KQSP+LRIN+LL+2).EQ.IHL(2).AND.
     +         IQ(KQSP+LRIN+LL+4).EQ.IRD)THEN
               CALL UCOPY2(IQ(KQSP+LRIN+LL+NLESS),
     +                     IQ(KQSP+LRIN+LL),NWL-LL+1)
               NWL=NWL-NLESS
               IQ(KQSP+LRIN+1)=IQ(KQSP+LRIN+1)-1
               IQ(KQSP+LRIN+2)=NWL
            ELSE
               LL=LL+NLESS
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
*
*           Delete list of allocated records
*
      IF(LFREE.NE.0)THEN
         CALL VZERO(IQ(KQSP+LFREE+1),IQ(KQSP+LFREE-1))
      ENDIF
*
*           Write back record 1
*
      CALL RZIODO(LUN,MAX(NWL,50),1,IQ(KQSP+LRIN+1),2)
      IQUEST(10)=IQ(KQSP+LRIN+1)
*
  99  RETURN
      END
