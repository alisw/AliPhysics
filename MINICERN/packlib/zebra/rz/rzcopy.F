*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:44  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZCOPY(CHPATH,KEYU,ICYCLE,KEYUN,CHOPT)
*
************************************************************************
*
*           Routine to copy an object from CHPATH or the whole tree to the CWD
* Input:
*   CHPATH  The pathname of the directory tree  which has to be copied to
*           the CWD
*   KEYU    KEY of the object to be copied from CHPATH
*   ICYCLE  Cycle number of the key to be copied
*   KEYUN   New value of the key in CWD (may be the same as KEYU)
*   CHOPT   Character string to specify various options
*   default ' ' copy the object with (KEYU,ICYCLE) from CHPATH to the CWD
*               If KEYUN already exists, a new cycle is created
*           'C' copy all cycles for the specified key
*           'K' copy all keys  (If 'C' option is given, copy all cycles)
*           'T' copy the complete tree CHPATH
*               When the option 'T' is given, by default only the highest
*               cycle of each key is copied
*               To copy all cycles use 'TC' option
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 07.05.86
*  Last mod: 14.05.92 Add CHOPT on call to RZFDIR
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzcycle.inc"
#if defined(CERNLIB_QMVAX)
#include "zebra/rzclun.inc"
#endif
      CHARACTER*(*) CHPATH,CHOPT
      DIMENSION KEYU(*),KEYUN(*)
      DIMENSION IOPTV(3),ISD(NLPATM),NSD(NLPATM),IHDIR(4)
      CHARACTER*16 CHFPAT(NLPATM)
      EQUIVALENCE (IOPTC,IOPTV(1)),(IOPTK,IOPTV(2))
     +           ,(IOPTT,IOPTV(3))
      LOGICAL COPY,RZSAME
*
*-----------------------------------------------------------------------
*
#if defined(CERNLIB_QMVAX)
#include "zebra/q_jbit.inc"
#endif
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      LOGLV=JBYT(IQ(KQSP+LTOP),15,3)-3
*
*          Save existing material (if any)
*
      CALL RZSAVE
*
      CALL UOPTC(CHOPT,'CKT',IOPTV)
*
*         Check if WRITE permission on file and directory
*
      IF(LQRS.EQ.0)GO TO 999
      IFLAG=0
      CALL RZMODS('RZCOPY',IFLAG)
      IF(IFLAG.NE.0)GO TO 999
*
*         Save CWD name
*
      CALL RZCDIR(CHWOLD,'R')
*
*         Load directory CHPATH
*
      CALL RZPATH(CHPATH)
      NLPAT0=NLPAT
      NLPAT1=NLPAT
      DO 1 I=1,NLPAT
   1  CHFPAT(I)=CHPAT(I)
      CALL RZFDIR('RZCOPY',LT,LFROM,' ')
      IF(LFROM.EQ.0)THEN
         IQUEST(1)=4
         GO TO 999
      ENDIF
      ISD(NLPAT1)=0
      NSD(NLPAT1)=IQ(KQSP+LFROM+KNSD)
      CALL SBIT0(IQ(KQSP+LFROM),IQDROP)
      LB    =IQ(KQSP+LT+KLB)
      LROLD =IQ(KQSP+LT+LB+1)
      LUNOLD=IQ(KQSP+LT-5)
      NKEYS=IQ(KQSP+LFROM+KNKEYS)
      NWKEY=IQ(KQSP+LFROM+KNWKEY)
*
*        Check if KEY descriptors matches
*
      IF(NWKEY.NE.IQ(KQSP+LCDIR+KNWKEY).OR.
     +   IQ(KQSP+LFROM+KKDES).NE.IQ(KQSP+LCDIR+KKDES))THEN
         IQUEST(1)=4
         IF(LOGLV.GE.-2) WRITE(IQLOG,1000)
 1000    FORMAT(' RZCOPY. Key descriptors do not match')
         GO TO 900
      ENDIF
*
      IF(IOPTT.NE.0)THEN
         LBANK=LCDIR
   5     IF(LBANK.NE.LTOP)THEN
            LBANK=LQ(KQSP+LBANK+1)
            IF(LBANK.EQ.LFROM)THEN
               IF(LOGLV.GE.-2) WRITE(IQLOG,3000)
 3000          FORMAT(' RZCOPY. Cannot copy mother tree in daughter')
               IQUEST(1)=4
               GO TO 900
            ENDIF
            GO TO 5
         ENDIF
      ENDIF
*
      IF(NKEYS.EQ.0)THEN
         IF(IOPTK.EQ.0.AND.IOPTT.EQ.0)GO TO 999
         GO TO 100
      ENDIF
*
*        Convert KEYU,KEYUN (If only one key to be copied)
*
      IF(IOPTK.EQ.0.AND.IOPTT.EQ.0)THEN
         DO 10 I=1,NWKEY
            IKDES=(I-1)/10
            IKBIT1=3*I-30*IKDES-2
            IF(JBYT(IQ(KQSP+LFROM+KKDES+IKDES),IKBIT1,3).LT.3)THEN
               KEY(I)=KEYU(I)
               KEY2(I)=KEYUN(I)
            ELSE
               CALL ZHTOI(KEYU(I),KEY(I),1)
               CALL ZHTOI(KEYUN(I),KEY2(I),1)
            ENDIF
  10     CONTINUE
      ENDIF
  15  IF(IOPTT.NE.0)THEN
         ISD(NLPAT1)=0
         NSD(NLPAT1)=IQ(KQSP+LFROM+KNSD)
       ENDIF
*
*        Loop on all keys of level 0
*
      DO 80 I=1,NKEYS
         LK=IQ(KQSP+LFROM+KLK)
         LKC=LK+(NWKEY+1)*(I-1)
         IF(IOPTK.EQ.0.AND.IOPTT.EQ.0)THEN
            DO 20 K=1,NWKEY
               IF(IQ(KQSP+LFROM+LKC+K).NE.KEY(K))GO TO 80
  20        CONTINUE
         ELSE
            DO 25 K=1,NWKEY
               KEY2(K)=IQ(KQSP+LFROM+LKC+K)
  25        CONTINUE
         ENDIF
         LCYC  =IQ(KQSP+LFROM+LKC)
         IF (KVSCYC.NE.0) THEN
*           IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.I) THEN
            IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.IQ(KQSP+LFROM+LKC+1)) THEN
               IQUEST(1) = 11
               GO TO 900
            ENDIF
         ENDIF
*
*        Store cycles in reverse order for 'C' option
*
         IF(IOPTC.NE.0)THEN
            IF(LCORD.EQ.0)THEN
               CALL MZBOOK(JQPDVS,LCORD,LTOP,-4,'RZCO',0,0,50,2,-1)
            ENDIF
            IQ(KQSP+LCORD+1)=0
  30        NORD=IQ(KQSP+LCORD+1)+1
            IF (KVSCYC.NE.0) THEN
               LCOLD = JBYT(IQ(KQSP+LFROM+LCYC+KPPCYC), 1,16)
            ELSE
               LCOLD = IQ(KQSP+LFROM+LCYC+KPPCYC)
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
               LCYC=IQ(KQSP+LCORD+IC+1)
               CALL RZCOP1(LUNOLD,LROLD,KEY2,IQ(KQSP+LFROM+LCYC+1))
               IF(IQUEST(1).NE.0) GO TO 900
  40        CONTINUE
         ELSE
  50        IF (KVSCYC.NE.0) THEN
               LCOLD = JBYT(IQ(KQSP+LFROM+LCYC+KPPCYC), 1,16)
            ELSE
               LCOLD = IQ(KQSP+LFROM+LCYC+KPPCYC)
            ENDIF
            ICY   = JBYT(IQ(KQSP+LFROM+LCYC+KCNCYC),21,12)
            COPY=ICYCLE.GE.ICY.OR.(ICYCLE.LE.0.AND.LCOLD.EQ.0).OR.
     +           IOPTT.NE.0.OR.IOPTK.NE.0
            IF(COPY)THEN
               CALL RZCOP1(LUNOLD,LROLD,KEY2,IQ(KQSP+LFROM+LCYC+1))
               IF(IQUEST(1).NE.0) GO TO 900
            ELSE
               IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
                 LCYC=LCOLD
                 GO TO 50
               ENDIF
            ENDIF
         ENDIF
         IF(IOPTK.EQ.0.AND.IOPTT.EQ.0)GO TO 900
*
  80  CONTINUE
*
 100  IF(IOPTT.EQ.0)GO TO 900
*
*             Copy subdirectories
*
 110  ISD(NLPAT1)=ISD(NLPAT1)+1
      IF(ISD(NLPAT1).LE.NSD(NLPAT1))THEN
         NLPAT1=NLPAT1+1
         LSF=IQ(KQSP+LFROM+KLS)
         IH=LSF+7*(ISD(NLPAT1-1)-1)
         CALL ZITOH(IQ(KQSP+LFROM+IH),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHFPAT(NLPAT1),16)
         DO 120 I=1,NLPAT1
 120     CHPAT(I)=CHFPAT(I)
         NLPAT=NLPAT1
         CALL RZFDIR('RZCOPY',LT,LFROM,' ')
         IF(LFROM.EQ.0)THEN
            IQUEST(1)=4
            GO TO 900
         ENDIF
         NKEYS=IQ(KQSP+LFROM+KNKEYS)
         NWKEY=IQ(KQSP+LFROM+KNWKEY)
         KTAGS=KKDES+(NWKEY-1)/10+1
         DO 130 I=2,NLPAT1
 130     CHCDIR(I)=CHFPAT(I)
         CALL RZPAFF(CHCDIR,NLPAT1-1,CHL)
         CALL RZCDIR(CHL,' ')
         CALL RZMDIR(CHFPAT(NLPAT1),NWKEY,'?',' ')
         IF(IQUEST(1).NE.0)GO TO 900
         CALL RZPAFF(CHCDIR,NLPAT1,CHL)
         CALL RZCDIR(CHL,' ')
         IF(IQ(KQSP+LCDIR-1).LT.2*NWKEY+KTAGS+20)THEN
            CALL RZEXPD('RZCOPY',100)
            IF(IQUEST(1).NE.0)GO TO 900
         ENDIF
         CALL UCOPY(IQ(KQSP+LFROM+KKDES),IQ(KQSP+LCDIR+KKDES),
     +              2*NWKEY+KTAGS-KKDES)
         CALL SBIT1(IQ(KQSP+LTOP),2)
         CALL SBIT1(IQ(KQSP+LCDIR),2)
         GO TO 15
      ELSE
         NLPAT1=NLPAT1-1
         IF(NLPAT1.GE.NLPAT0)THEN
            LUP=LQ(KQSP+LFROM+1)
            CALL MZDROP(JQPDVS,LFROM,' ')
            LFROM=LUP
            GO TO 110
         ENDIF
      ENDIF
*
 900  IRCOD = IQUEST(1)
      IF(LCORD.NE.0)THEN
         CALL MZDROP(JQPDVS,LCORD,' ')
         LCORD=0
      ENDIF
      IF(LRIN.NE.0)THEN
         CALL MZDROP(JQPDVS,LRIN ,' ')
         LRIN=0
      ENDIF
      CALL RZCDIR(CHWOLD,' ')
      IF(LFROM.NE.LCDIR)CALL SBIT1(IQ(KQSP+LFROM),IQDROP)
      IQUEST(1) = IRCOD
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
      IF(LUNOLD.NE.LUN)THEN
         IF(JBIT(IQ(KQSP+LFROM),4).NE.0)UNLOCK(UNIT=LUNOLD)
      ENDIF
#endif
*
 999  RETURN
      END
