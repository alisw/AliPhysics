*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:10  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZTOF1(LUNFZ,IOPTC)
*
************************************************************************
*
*        Copy current directory to a sequential FZ file
* Input:
*   LUNFZ   Logical unit number of the FZ sequential access file
*   IOPTC   0 save only the highest cycle to LUNFZ
*           1 save all cycles
*
* Called by <RZTOFZ>
*
*  Author  : R.Brun DD/US/PD
*  Written : 14.05.86
*  Last mod: 04.10.90
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
*            Fill header for directory
*
      NKEYS=IQ(KQSP+LCDIR+KNKEYS)
      NWKEY=IQ(KQSP+LCDIR+KNWKEY)
      LB   = IQ(KQSP+LTOP+KLB)
      LREC = IQ(KQSP+LTOP+LB+1)
      IDECK=0
      IF(NWKEY.EQ.2.AND.LREC.EQ.128)THEN
         KTAGS=KKDES+1
         CALL ZITOH(IQ(KQSP+LCDIR+KTAGS),KEY,2)
         CALL UCTOH('DECKNAME',KEY2,4,8)
         IF(KEY(1).EQ.KEY2(1).AND.KEY(2).EQ.KEY2(2))THEN
            IDECK=1
         ENDIF
      ENDIF
      NH   =NWKEY+3
      IHEAD(1)=1
      IHEAD(2)=NLPAT
      IHEAD(3)=0
      CALL UCOPY(IQ(KQSP+LCDIR+1),KEY,KNSIZE-3)
*
*            Write directory header
*
      CALL FZOUT(LUNFZ,JQPDVS,0,1,'Z',1,KNSIZE,IHEAD)
      IF(IQUEST(1).NE.0)GO TO 99
*
*           Loop on all keys of level 0
*
      IHEAD(1)=0
      DO 80 I=1,NKEYS
         LK=IQ(KQSP+LCDIR+KLK)
         LKC=LK+(NWKEY+1)*(I-1)
         DO 25 K=1,NWKEY
            KEY(K)=IQ(KQSP+LCDIR+LKC+K)
  25     CONTINUE
         LCYC  =IQ(KQSP+LCDIR+LKC)
         IF (KVSCYC.NE.0) THEN
*           IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.I) THEN
            IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.IQ(KQSP+LCDIR+LKC+1)) THEN
               IQUEST(1) = 11
               GO TO 99
            ENDIF
         ENDIF
*
*               Store cycles in reverse order for 'C' option
*
         IF(IOPTC.NE.0)THEN
            IF(LCORD.EQ.0)THEN
               CALL MZBOOK(JQPDVS,LCORD,LTOP,-4,'RZCO',0,0,50,2,-1)
            ENDIF
            IQ(KQSP+LCORD+1)=0
  30        NORD=IQ(KQSP+LCORD+1)+1
            IF (KVSCYC.EQ.0) THEN
               LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC+KPPCYC), 1,16)
            ELSE
               LCOLD = IQ(KQSP+LCDIR+LCYC+KPPCYC)
            ENDIF
            IF(NORD.GT.IQ(KQSP+LCORD-1))THEN
               CALL MZPUSH(JQPDVS,LCORD,0,50,'I')
            ENDIF
            IQ(KQSP+LCORD+1)=NORD
            IQ(KQSP+LCORD+NORD+1)=LCYC
            IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
               LCYC=LCOLD
               GO TO 30
            ENDIF
            DO 40 IC=NORD,1,-1
               LCYC= IQ(KQSP+LCORD+IC+1)
               ICY = JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
               CALL RZIN(JQPDVS,LRZ0,-1,I,ICY,'S')
               IF(IQUEST(1).NE.0)GO TO 99
               LFROM=LQ(KQSP+LRZ0-1)
               IF(IDECK.NE.0)THEN
                  NDATA=IQ(KQSP+LFROM-1)
                  NT=NDATA/5 +40
                  CALL MZBOOK(JQPDVS,LTEMP,LFROM,0,'TEMP',0,0,NT,1,-1)
                  NTOT=0
                  NEW=1
  35              IF(NTOT.LT.NDATA)THEN
                     IF(NEW.GT.NT-1)THEN
                        CALL MZPUSH(JQPDVS,LTEMP,0,50,'I')
                        NT=NT+50
                     ENDIF
                     CALL RZLIND(IQ(KQSP+LFROM+1),NTOT,
     +                           IQ(KQSP+LTEMP+1),NEW)
                     NEW=NEW+1
                     GO TO 35
                  ENDIF
               ENDIF
               IHEAD(2)=IC
               IHEAD(3)=IQ(KQSP+LCDIR+LCYC+KFLCYC)
               IQUEST(1)=0
               CALL FZOUT(LUNFZ,JQPDVS,LFROM,1,'L',2,NH,IHEAD)
               IF(IQUEST(1).NE.0)GO TO 90
               CALL MZDROP(JQPDVS,LFROM,'L')
               LFROM=0
  40        CONTINUE
         ELSE
  50        ICY   =JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
            CALL RZIN(JQPDVS,LRZ0,-1,I,ICY,'S')
            IF(IQUEST(1).NE.0)GO TO 99
            LFROM=LQ(KQSP+LRZ0-1)
            IF(IDECK.NE.0)THEN
               NDATA=IQ(KQSP+LFROM-1)
               NT=NDATA/5 +40
               CALL MZBOOK(JQPDVS,LTEMP,LFROM,0,'TEMP',0,0,NT,1,-1)
               NTOT=0
               NEW=1
  55           IF(NTOT.LT.NDATA)THEN
                  IF(NEW.GT.NT-1)THEN
                     CALL MZPUSH(JQPDVS,LTEMP,0,50,'I')
                     NT=NT+50
                  ENDIF
                  CALL RZLIND(IQ(KQSP+LFROM+1),NTOT,
     +                        IQ(KQSP+LTEMP+1),NEW)
                  NEW=NEW+1
                  GO TO 55
               ENDIF
            ENDIF
            IHEAD(2)=1
            IHEAD(3)=IQ(KQSP+LCDIR+LCYC+KFLCYC)
            IQUEST(1)=0
            CALL FZOUT(LUNFZ,JQPDVS,LFROM,1,'L',2,NH,IHEAD)
            IF(IQUEST(1).NE.0)GO TO 90
            CALL MZDROP(JQPDVS,LFROM,'L')
            LFROM=0
         ENDIF
*
  80  CONTINUE
*
*           Write directory trailer
*
      CALL FZOUT(LUNFZ,JQPDVS,0,0,'Z',1,1,77)
*
  90  IF(LFROM.GT.0)THEN
         IQ1=IQUEST(1)
         CALL MZDROP(JQPDVS,LFROM,'L')
         IQUEST(1)=IQ1
      ENDIF
*
  99  RETURN
      END
