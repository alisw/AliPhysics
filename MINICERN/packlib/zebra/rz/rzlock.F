*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:01  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLOCK(CHLOCK)
*
************************************************************************
*
*        To lock the CWD
* Input:
*   CHLOCK  Character variable  identifying the owner  of the  lock (e.g.
*           specifying the name of the user, his computer identifier,...)
*           This parameter is used to avoid two users,  who have both the
*           write password for  a directory,  trying to change  it at the
*           same time.    CHLOCK is also useful  in the case of  a system
*           crash while a directory was locked.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.05.86
*  Last mod: 04.10.90
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
      CHARACTER*(*) CHLOCK
      DIMENSION IHL(2),KHL(2)
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      IF(LTOP.EQ.0)GO TO 99
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
*** 9010    FORMAT(' RZLOCK. No authorisation to write in that directory')
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
         IF(NTRY.LT.100.AND.IQUEST(1).EQ.0)GO TO 10
         IF(LOGLV.GE.-2) WRITE(IQLOG,1000)
 1000    FORMAT(' RZLOCK. Cannot lock that directory')
         IQUEST(1)=1
         GO TO 99
      ENDIF
      IQ(KQSP+LRIN+3)=1
      CALL RZIODO(LUN,MAX(NWL,50),1,IQ(KQSP+LRIN+1),2)
*
      IQ(KQSP+LRIN+3)=0
      NLOCK=IQ(KQSP+LRIN+1)
      NREC=IQ(KQSP+LTOP+KQUOTA)
      LDC=IQ(KQSP+LCDIR+KLD)
      IRD=IQ(KQSP+LCDIR+LDC+1)
*
*           Read fresh Top directory
*
      IF(LTOP.NE.LCDIR)CALL RZRTOP
*
*          Check that directory is not already locked
*
      IF(NLOCK.GT.0)THEN
         LL=4
  20     IF(IQ(KQSP+LRIN+LL).NE.0)THEN
*
*          Check mother directories
*
            IMOT=0
            IRCUR=IQ(KQSP+LRIN+LL+4)
            IF(IRCUR.EQ.2.OR.IRCUR.EQ.IRD)IMOT=1
            IF(IMOT.NE.0)THEN
               CALL ZITOH(IQ(KQSP+LRIN+LL+1),KHL,2)
               IF(LOGLV.GE.-2) WRITE(IQLOG,3000)KHL
 3000          FORMAT(' RZLOCK. Directory already locked by ',2A4)
               IQUEST(1)=2
               GO TO 90
            ELSE
               LL=LL+IQ(KQSP+LRIN+LL)
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
*
*          Fill 'free' bank with allocated records as a function
*          of quota
*
      NFREE=0
      NRUSED=IQ(KQSP+LCDIR+KRUSED)
      NDATA=IQ(KQSP+LFREE-1)
      LB=IQ(KQSP+LTOP+KLB)
      IFR=2
      DO 50 I=3,NREC
         LL=4
  30     NMORE=IQ(KQSP+LRIN+LL)
         IF(NMORE.NE.0)THEN
            ND=IQ(KQSP+LRIN+LL+5)
            DO 40 J=1,ND
               IR1=IQ(KQSP+LRIN+LL+2*J+4)
               IRL=IQ(KQSP+LRIN+LL+2*J+5)
               IF(I.GE.IR1.AND.I.LE.IRL)GO TO 50
  40        CONTINUE
            LL=LL+NMORE
            GO TO 30
         ENDIF
*
         IWORD=(I-1)/32+1
         IBIT=I-32*(IWORD-1)
         IF(JBIT(IQ(KQSP+LTOP+LB+2+IWORD),IBIT).EQ.0)THEN
            NRUSED=NRUSED+1
            IF(NRUSED.GT.IQ(KQSP+LCDIR+KQUOTA))GO TO 60
            IF(IQ(KQSP+LFREE+IFR).EQ.0)THEN
               NFREE=NFREE+1
               IQ(KQSP+LFREE+1)=NFREE
               IQ(KQSP+LFREE+IFR)=I
               IQ(KQSP+LFREE+IFR+1)=I
            ELSE
               IF(I.EQ.IQ(KQSP+LFREE+IFR+1)+1)THEN
                  NFREE=IQ(KQSP+LFREE+1)
                  IF(NFREE.EQ.0)NFREE=1
                  IQ(KQSP+LFREE+IFR+1)=I
               ELSE
                  NFREE=NFREE+1
                  IQ(KQSP+LFREE+1)=NFREE
                  IF(2*NFREE+3.GT.NDATA)THEN
                     CALL MZPUSH(JQPDVS,LFREE,0,20,'I')
                     NDATA=NDATA+20
                  ENDIF
                  IFR=IFR+2
                  IQ(KQSP+LFREE+IFR)=I
                  IQ(KQSP+LFREE+IFR+1)=I
               ENDIF
            ENDIF
         ENDIF
  50  CONTINUE
*
*          Build new lock
*
  60  NMORE=2*NFREE+6
      IF(NFREE.LE.0)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,3100)
 3100    FORMAT(' RZLOCK. Cannot allocate free records -',
     +          ' RZ quota for this file has been reached.')
         IQUEST(1)=3
         GO TO 90
      ENDIF
      IF(NWL.EQ.0)NWL=4
      IF(NWL+NMORE.GT.LREC)THEN
         NF=NFREE
         NFREE=(LREC-NWL-6)/2
         IQ(KQSP+LFREE+1)=NFREE
         IF(NFREE.GT.0)THEN
            CALL UCOPY(IQ(KQSP+LFREE+2*NF),IQ(KQSP+LFREE+2*NFREE),2)
            IQ(KQSP+LFREE+2*NFREE+2)=0
            NMORE=2*NFREE+6
            IF(LOGLV.GE.-2) WRITE(IQLOG,4000)
 4000       FORMAT(' RZLOCK. Cannot allocate all free records')
         ELSE
            IF(LOGLV.GE.-2) WRITE(IQLOG,4100)
 4100       FORMAT(' RZLOCK. Data base is too fragmented')
            IQUEST(1)=1
            IQ(KQSP+LFREE+1)=0
            GO TO 90
         ENDIF
      ENDIF
      IQ(KQSP+LRIN+NWL)=NMORE
      IQ(KQSP+LRIN+NWL+1)=IHL(1)
      IQ(KQSP+LRIN+NWL+2)=IHL(2)
      IQ(KQSP+LRIN+NWL+3)=0
      CALL RZDATE(IQ(KQSP+LRIN+NWL+3),IDATE,ITIME,2)
      IQ(KQSP+LRIN+NWL+4)=IQ(KQSP+LCDIR+LDC+1)
      CALL UCOPY(IQ(KQSP+LFREE+1),IQ(KQSP+LRIN+NWL+5),2*NFREE+1)
      NWL=NWL+NMORE
      IQ(KQSP+LRIN+NWL)=0
      IQ(KQSP+LRIN+1)=IQ(KQSP+LRIN+1)+1
      IQ(KQSP+LRIN+2)=NWL
*
*          Reset the lock and write record 1
*
  90  CALL RZIODO(LUN,MAX(NWL,50),1,IQ(KQSP+LRIN+1),2)
      IQUEST(10)=IQ(KQSP+LRIN+1)
*
  99  RETURN
      END
