*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:00  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLLOK
*
************************************************************************
*
*       Routine to print current active locks
*
*
* Called by RZFILE
*
*  Author  : R.Brun DD/US/PD
*  Written : 08.09.89
*  Last mod: 08.09.89
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      DIMENSION IDIR(5,10),KHL(2)
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      IF(LTOP.EQ.0)GO TO 99
*
*           Read locking record
*
      IF(LRIN.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LRIN,LTOP,-7,'RZIN',0,0,LREC+1,2,-1)
         IQ(KQSP+LRIN-5)=IQ(KQSP+LTOP-5)
      ENDIF
      NWL =50
      NTRY=0
  10  CALL RZIODO(LUN,NWL,1,IQ(KQSP+LRIN+1),1)
      IF(IQUEST(1).NE.0)GO TO 90
      IF(IQ(KQSP+LRIN+2).GT.NWL)THEN
         NWL=IQ(KQSP+LRIN+2)
         GO TO 10
      ENDIF
      IQ(KQSP+LTOP+KIRIN)=0
      IF(IQ(KQSP+LRIN+3).NE.0)THEN
         NWL=50
         NTRY=NTRY+1
#if defined(CERNLIB_QMVAX)
         IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
         CALL LIB$WAIT(0.1)
#endif
         IF(NTRY.LT.100)GO TO 10
         IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,1000)
 1000    FORMAT(' RZLLOK. Cannot get locking record')
         IQUEST(1)=1
         GO TO 90
      ENDIF
*
      NLOCK=IQ(KQSP+LRIN+1)
      IF(NLOCK.LE.0)GO TO 99
      LL=4
  20  IF(IQ(KQSP+LRIN+LL).NE.0)THEN
         IRD=IQ(KQSP+LRIN+LL+4)
         CALL ZITOH(IQ(KQSP+LRIN+LL+1),KHL,2)
         IDTIME=IQ(KQSP+LRIN+LL+3)
         CALL RZDATE(IDTIME,IDATE,ITIME,1)
         NLEVEL=11
  30     NLEVEL=NLEVEL-1
         CALL RZIODO(LUN,5,IRD,IDIR(1,NLEVEL),1)
         IF(IQUEST(1).NE.0)GO TO 90
         CALL ZITOH(IDIR(1,NLEVEL),IDIR(1,NLEVEL),4)
         IRD=IDIR(5,NLEVEL)
         IF(IRD.GT.0)GO TO 30
*
         NL=11-NLEVEL
         CALL UCOPY2(IDIR(1,NLEVEL),IDIR(1,1),NL*5)
         DO 40 I=1,NL
            CALL UHTOC(IDIR(1,I),4,CHPAT(I),16)
  40     CONTINUE
         CALL RZPAFF(CHPAT,NL,CHL)
         WRITE(IQPRNT,2000)KHL,IDATE,ITIME,CHL(1:70)
 2000    FORMAT(' LOCK-ID < ',2A4,'> on ',I6,'/',I4,' for directory ',A)
*
         LL=LL+IQ(KQSP+LRIN+LL)
         GO TO 20
      ENDIF
  90  CONTINUE
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
*
  99  RETURN
      END
