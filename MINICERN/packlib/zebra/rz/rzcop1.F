*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:43  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZCOP1(LUNOLD,LROLD,KEY,IOLD)
*
************************************************************************
*
*          Copy one (KEY,CYCLE) from LFROM to the CWD
* Input:
*   LUNOLD  Logical unit number of the file from which the copy is made
*   LROLD   Record size of LUNOLD
*   KEY     Identification (array) of the key to created in the CWD
*   IOLD    Array of 4 words describing the cycle which is copied
*
* Called by RZCOPY
*
*  Author  : R.Brun DD/US/PD
*  Written : 07.05.86
*  Last mod: 01.09.92 Dave Morrison (MIT) handle append mode bit
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
      DIMENSION KEY(*),IOLD(4)
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

*
*         Get last record written in that directory
*         Create buffer bank
*
      IF (KVSCYC.EQ.0) THEN
         IR1OLD = JBYT(IOLD(KFRCYC),17,16)
         IP1OLD = JBYT(IOLD(KORCYC), 1,16)
         IR2OLD = JBYT(IOLD(KSRCYC),17,16)
         NDATA  = JBYT(IOLD(KNWCYC), 1,20)
         IFORM  = JBYT(IOLD(KFLCYC), 1, 3)
      ELSE
         IR1OLD = IOLD(KFRCYC)
         IP1OLD = JBYT(IOLD(KORCYC), 1,20)
         IR2OLD = IOLD(KSRCYC)
         NDATA  = IOLD(KNWCYC)
         IFORM  = JBYT(IOLD(KFLCYC), 1, 3)
      ENDIF
      LROUT  = LQ(KQSP+LTOP-6)
      IROUT  = IQ(KQSP+LTOP+KIROUT)
      IRLOUT = IQ(KQSP+LCDIR+KRLOUT)
      IP1    = IQ(KQSP+LCDIR+KIP1)
      IF(LROUT.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LROUT,LTOP,-6,'RZOU',0,0,LREC+1,2,-1)
         IQ(KQSP+LROUT-5)=LUN
         IROUT=0
         IP1=1
      ENDIF
      IF(IROUT.NE.IRLOUT.AND.IRLOUT.NE.0)THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),1)
         IF(IQUEST(1).NE.0)GO TO 999
#if defined(CERNLIB_QMVAX)
         IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
         IROUT=IRLOUT
         IQ(KQSP+LTOP+KIROUT)=IROUT
         IP1=IQ(KQSP+LCDIR+KIP1)
         IF(IQ(KQSP+LTOP+KIRIN).EQ.IROUT)IQ(KQSP+LTOP+KIRIN) = 0
      ENDIF
*
*         Is directory big enough to accomodate new cycle ?
*
      IF(IQ(KQSP+LCDIR+KNFREE).LT.IQ(KQSP+LCDIR+KNWKEY)+4*KLCYCL+1)THEN
         CALL RZEXPD('RZCOPY',10*(IQ(KQSP+LCDIR+KNWKEY)+KLCYCL+1))
         IF(IQUEST(1).NE.0) GO TO 999
      ENDIF
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
      LC = IQ(KQSP+LCDIR+KLC)
      NWFREE=IQ(KQSP+LCDIR+KNFREE)
*
*         Compute how many records
*         are necessary to write data structure.
*
      NLEFT=LREC-IP1+1
      IF(NDATA.LE.NLEFT)THEN
         N1=NDATA
         NR=0
      ELSE
         N1=NLEFT
         NR=(NDATA-NLEFT-1)/LREC + 1
      ENDIF
      IF(IRLOUT.EQ.0)NR=NR+1
      IF(NR.GT.0)THEN
         CALL RZALLO('RZCOPY',NR,IALLOC)
         IF(IALLOC.EQ.0) GO TO 999
         IF(IRLOUT.EQ.0)IRLOUT=IALLOC
      ENDIF
*
*         Search if KEY is already entered
*
      NKEYS = IQ(KQSP+LCDIR+KNKEYS)
      NWKEY = IQ(KQSP+LCDIR+KNWKEY)
      IQUEST(7)=NKEYS
      IQUEST(8)=NWKEY
*
      IF(NKEYS.GT.0)THEN
         DO 20 I=1,NKEYS
            LKC=LK+(NWKEY+1)*(I-1)
            DO 10 K=1,NWKEY
               IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GO TO 20
  10        CONTINUE
            LCOLD  = IQ(KQSP+LCDIR+LKC)
            IF (KVSCYC.NE.0) THEN
*              IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.I) THEN
               IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.
     +             IQ(KQSP+LCDIR+LKC+1)) THEN
                  IQUEST(1) = 11
                  GO TO 999
               ENDIF
            ENDIF
            ICOLD  = JBYT(IQ(KQSP+LCDIR+LCOLD+KCNCYC),21,12)
            ICYCLE = ICOLD+1
*           IKYV   = I
            IKYV   = IQ(KQSP+LCDIR+LKC+1)
            GO TO 50
  20     CONTINUE
      ENDIF
*
*         New KEY, append to the list
*
      NWFREE=NWFREE-NWKEY-1
      IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)+1
      LCOLD=0
      LKC=LF
      LF=LF+NWKEY+1
      DO 30 I=1,NWKEY
         IQ(KQSP+LCDIR+LKC+I)=KEY(I)
  30  CONTINUE
      ICYCLE = 1
*     IKYV   = IQ(KQSP+LCDIR+KNKEYS)
      IKYV   = IQ(KQSP+LCDIR+LKC+1)
*
*         Create a new cycle
*
  50  LKCSV = IQ(KQSP+LCDIR+LKC)
      LC    = LC-KLCYCL
      NWFREE= NWFREE-KLCYCL
      IQ(KQSP+LCDIR+LKC) = LC
      IQ(KQSP+LCDIR+LC+KPPCYC) = LCOLD
      IQ(KQSP+LCDIR+LC+KFLCYC)=0
      CALL RZDATE (IQ(KQSP+LCDIR+LC+KFLCYC),IDATE,ITIME,2)
c
c     DPM: Added this to handle append mode.
c
      IBIT4 = JBIT(IOLD(KFLCYC),4)
      IF(IBIT4.EQ.1)CALL SBIT1(IQ(KQSP+LCDIR+LC+KFLCYC),4)
c
c     DPM: End of changes.
c
      CALL SBYT (IFORM,IQ(KQSP+LCDIR+LC+KFLCYC),1,3)
      IQ(KQSP+LCDIR+LC+KORCYC) = IP1
      IQ(KQSP+LCDIR+LC+KNWCYC) = NDATA
      CALL SBYT(ICYCLE,IQ(KQSP+LCDIR+LC+KCNCYC),21,12)
      IF (KVSCYC.EQ.0) THEN
         IF(N1.LT.NDATA)CALL SBYT(IALLOC,IQ(KQSP+LCDIR+LC+KSRCYC),17,16)
         CALL SBYT(IRLOUT,IQ(KQSP+LCDIR+LC+KFRCYC),17,16)
      ELSE
         IF (N1.LT.NDATA) THEN
            IQ(KQSP+LCDIR+LC+KSRCYC) = IALLOC
         ELSE
            IQ(KQSP+LCDIR+LC+KSRCYC) = 0
         ENDIF
         IQ(KQSP+LCDIR+LC+KFRCYC) = IRLOUT
         IQ(KQSP+LCDIR+LC+KKYCYC) = IKYV
      ENDIF
      IQUEST(3)=IRLOUT
      IQUEST(4)=IP1
      IQUEST(5)=0
      IQUEST(6)=ICYCLE
      IQUEST(11)=NDATA
*
*          Copy records
*           Start filling current block
*
      IF(LRIN.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LRIN,LTOP,-7,'RZIN',0,0,LROLD+1,2,-1)
         IQ(KQSP+LRIN-5)=IQ(KQSP+LTOP-5)
      ELSE
         NOLD=IQ(KQSP+LRIN-1)
         IF(NOLD.LT.LROLD)THEN
            CALL MZPUSH(JQPDVS,LRIN,0,LROLD-NOLD,'I')
         ENDIF
      ENDIF
      CALL RZIODO(LUNOLD,LROLD,IR1OLD,IQ(KQSP+LRIN+1),1)
      IF(IQUEST(1).NE.0) GO TO 900
      NWC=N1
      IRN=IR2OLD-1
  55  IF(NWC.GT.LROLD-IP1OLD+1)THEN
         NWC=LROLD-IP1OLD+1
         CALL UCOPY(IQ(KQSP+LRIN+IP1OLD),IQ(KQSP+LROUT+IP1),NWC)
         IP1OLD=1
         IF(NWC.LT.N1)THEN
            IRN=IRN+1
            CALL RZIODO(LUNOLD,LROLD,IRN,IQ(KQSP+LRIN+1),1)
            IF(IQUEST(1).NE.0) GO TO 900
            NWC=N1-NWC
            GO TO 55
         ENDIF
      ELSE
         CALL UCOPY(IQ(KQSP+LRIN+IP1OLD),IQ(KQSP+LROUT+IP1),NWC)
         IP1OLD=IP1OLD+NWC
      ENDIF
      IF(IP1.EQ.1)THEN
         IRLOUT=IALLOC
         IROUT=IRLOUT
      ENDIF
      IP1=IP1+N1
      IF(IP1.GT.LREC)THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),2)
         IF(IQUEST(1).NE.0) GO TO 900
         IF(IP1.GT.NDATA)IRLOUT=0
         IP1=1
      ENDIF
      IF(N1.LT.NDATA)THEN
         IQUEST(5)=IALLOC
         IQUEST(2)=NR+1
         DO 60 I=1,NR
            IP1=1
            NW=NDATA-N1
            IF(NW.GT.LREC)NW=LREC
            NWC=NW
  57        IF(NWC.GT.LROLD-IP1OLD+1)THEN
               NWC=LROLD-IP1OLD+1
               CALL UCOPY(IQ(KQSP+LRIN+IP1OLD),IQ(KQSP+LROUT+IP1),NWC)
               IP1OLD=1
               IF(NWC.LT.N1)THEN
                  IRN=IRN+1
                  CALL RZIODO(LUNOLD,LROLD,IRN,IQ(KQSP+LRIN+1),1)
                  IF(IQUEST(1).NE.0) GO TO 900
                  NWC=NW-NWC
                  GO TO 57
               ENDIF
            ELSE
               CALL UCOPY(IQ(KQSP+LRIN+IP1OLD),IQ(KQSP+LROUT+IP1),NWC)
               IP1OLD=IP1OLD+NWC
            ENDIF
            IF(NW.EQ.LREC)THEN
               CALL RZIODO(LUN,LREC,IROUT,IQ(KQSP+LROUT+1),2)
               IF(IQUEST(1).NE.0) GO TO 900
            ELSE
               IRLOUT=IALLOC+I-1
               IROUT=IRLOUT
            ENDIF
            IP1=IP1+NW
            N1=N1+NW
  60     CONTINUE
      ENDIF
*
*         Update internal pointers in the directory
*
      IQ(KQSP+LTOP+KIROUT)=IROUT
      IQUEST(9)=IQ(KQSP+LCDIR+KQUOTA)-IQ(KQSP+LCDIR+KRUSED)
      IQ(KQSP+LCDIR+KRUSED)=IQ(KQSP+LCDIR+KRUSED)+NR
      NWUSED=IQ(KQSP+LCDIR+KWUSED)+NDATA
      IF(NWUSED.GT.1000000)THEN
         IQ(KQSP+LCDIR+KMEGA)=IQ(KQSP+LCDIR+KMEGA)+1
         IQ(KQSP+LCDIR+KWUSED)=NWUSED-1000000
      ELSE
         IQ(KQSP+LCDIR+KWUSED)=NWUSED
      ENDIF
      IQ(KQSP+LCDIR+KRLOUT)=IRLOUT
      IQ(KQSP+LCDIR+KIP1)=IP1
      IQ(KQSP+LCDIR+KNFREE)=NWFREE
      IQ(KQSP+LCDIR+KLF)=LF
      IQ(KQSP+LCDIR+KLC)=LC
*
*         Mark used records
*
      IF(NR.GT.0)THEN
         CALL RZUSED(NR,IALLOC)
      ENDIF
      GO TO 999
*         Reset internal pointers in case of I/O problem
*
 900  IF(ICYCLE.EQ.1)THEN
        IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)-1
      ENDIF
      IQ(KQSP+LCDIR+LKC)=LKCSV
*
 999  RETURN
      END
