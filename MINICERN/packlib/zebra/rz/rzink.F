*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:55  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZINK(KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*         To find and decode KEYU,ICYCLE
* Input:
*   KEYU    Keyword vector of the information to be read
*   ICYCLE  Cycle number of the key to be read
*           ICYCLE > highest cycle number means read the highest cycle
*           ICYCLE = 0 means read the lowest cycle
*   CHOPT   Character variable specifying the options selected.
*           data structure
*             default
*                   Same as 'D' below
*             'A'   Read continuation of the previously read data structure
*                   with identifier KEYU,ICYCLE
*                   Given that option implies that the record was written with
*                   the same option by a call to RZOUT.
*             'C'   Provide   information   about   the   cycle   numbers
*                   associated with KEY.
*                   The  total number  of  cycles  and the  cycle  number
*                   identifiers of the 19 highest  cycles are returned in
*                   IQUEST(50) and IQUEST(51..89) respectively
*             'D'   Read the  Data structure  with the  (key,cycle)  pair
*                   specified.
*             'N'   Read the neighbouring. keys (i.e. those preceding and
*                   following KEY).
*                   The  key-vectors of  the previous  and  next key  are
*                   available   respectively   as   IQUEST(31..35)    and
*                   IQUEST(41..45), see below.
*             'R'   Read data into existing bank at LSUP,JBIAS
*             'S'   KEYU(1) contains the key serial number
*                   IQUEST(20)= serial number of the key in directory
*                   IQUEST(21..20+NWKEY)=KEY(1....NWKEY)
*
* Called by RZIN,RZVIN
*
*  Author  : R.Brun DD/US/PD
*  Written : 09.05.86
*  Last mod: 11.09.89
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 23.03.95 J.Shiers - check on K/C blocks is on KEY(1)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzcout.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*(*) CHOPT
      DIMENSION KEYU(*)
      EQUIVALENCE (IOPTA,IQUEST(91)), (IOPTC,IQUEST(92))
     +,    (IOPTD,IQUEST(93)), (IOPTN,IQUEST(94)), (IOPTR,IQUEST(95))
     +,    (IOPTS,IQUEST(96))
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
      IQUEST(1)=0
      CALL UOPTC(CHOPT,'ACDNRS',IQUEST(91))
*
*           Search KEY and CYCLE
*
      LK=IQ(KQSP+LCDIR+KLK)
      NKEYS=IQ(KQSP+LCDIR+KNKEYS)
      NWKEY=IQ(KQSP+LCDIR+KNWKEY)
      IQUEST(7)=NKEYS
      IQUEST(8)=NWKEY
      IF(NKEYS.EQ.0)GO TO 90
*
      IF(IOPTS.NE.0)THEN
         IK1=KEYU(1)
         IK2=IK1
         IF(IK1.GT.NKEYS.OR.IK1.LE.0)THEN
            IQUEST(1)=1
            IQUEST(2)=IK1
            RETURN
         ENDIF
      ELSE
         IK1=1
         IK2=NKEYS
         DO 5 I=1,NWKEY
            IKDES=(I-1)/10
            IKBIT1=3*I-30*IKDES-2
            IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
               KEY(I)=KEYU(I)
            ELSE
               CALL ZHTOI(KEYU(I),KEY(I),1)
            ENDIF
   5     CONTINUE
      ENDIF
      DO 30 I=IK1,IK2
         LKC=LK+(NWKEY+1)*(I-1)
         IF(IOPTS.EQ.0)THEN
            DO 10 K=1,NWKEY
               IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GO TO 30
  10        CONTINUE
         ELSE
            DO 15 K=1,NWKEY
               IF(K.LT.10)THEN
                  IKDES=(K-1)/10
                  IKBIT1=3*K-30*IKDES-2
                  IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
                     IQUEST(20+K)=IQ(KQSP+LCDIR+LKC+K)
                  ELSE
                     CALL ZITOH(IQ(KQSP+LCDIR+LKC+K),IQUEST(20+K),1)
                  ENDIF
               ENDIF
  15        CONTINUE
         ENDIF
         IQUEST(20)=I
         LCYC=IQ(KQSP+LCDIR+LKC)
         IF (KVSCYC.NE.0) THEN
*           IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.I) THEN
*
*    Check should be on content of KEY(1)
*
            IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.IQ(KQSP+LCDIR+LKC+1)) THEN
               IQUEST(1) = 11
               GO TO 99
            ENDIF
         ENDIF
         NC=0
  20     NC=NC+1
         ICY = JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
         IF(ICY.EQ.ICYCLE)GO TO 50
         IF(NC.EQ.1.AND.ICYCLE.GT.ICY)GO TO 50
         IF (KVSCYC.EQ.0) THEN
            LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC+KPPCYC),1,16)
         ELSE
            LCOLD = IQ(KQSP+LCDIR+LCYC+KPPCYC)
         ENDIF
         IF(LCOLD.EQ.0.AND.LCOLD.NE.LCYC.AND.ICYCLE.EQ.0)GO TO 50
         LCYC=LCOLD
         IF(LCYC.NE.0)GO TO 20
         GO TO 90
  30  CONTINUE
      GO TO 90
*
*           Cycle has been found
*           Read record descriptor
*
  50  IF (KVSCYC.EQ.0) THEN
         IR1   = JBYT(IQ(KQSP+LCDIR+LCYC+KFRCYC),17,16)
         IR2   = JBYT(IQ(KQSP+LCDIR+LCYC+KSRCYC),17,16)
         IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,16)
         NW    = JBYT(IQ(KQSP+LCDIR+LCYC+KNWCYC), 1,20)
      ELSE
         IR1   = IQ(KQSP+LCDIR+LCYC+KFRCYC)
         IR2   = IQ(KQSP+LCDIR+LCYC+KSRCYC)
         IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,20)
         NW    = IQ(KQSP+LCDIR+LCYC+KNWCYC)
      ENDIF
      N1    = NW
      IQUEST(2)=1
      IF(IR2.NE.0)IQUEST(2)=(NW-N1-1)/LREC+2
      IQUEST(3)=IR1
      IQUEST(4)=IP1
      IQUEST(5)=IR2
      IQUEST(6)=ICY
      IQUEST(12)=NW
      IQUEST(14)=IQ(KQSP+LCDIR+LCYC+1)
      IQUEST(15)=LCYC
C
C           C option given
C
      IF(IOPTC.NE.0)THEN
         IQUEST(50)=0
         LC1=LCYC
  51     IQUEST(50)=IQUEST(50)+1
         IF (KVSCYC.EQ.0) THEN
            LCOLD = JBYT(IQ(KQSP+LCDIR+LC1+KPPCYC),1,16)
         ELSE
            LCOLD = IQ(KQSP+LCDIR+LC1+KPPCYC)
         ENDIF
         IF(IQUEST(50).LE.19)THEN
            NC=IQUEST(50)
            IQUEST(50+NC)=JBYT(IQ(KQSP+LCDIR+LC1+KCNCYC),21,12)
            IQUEST(70+NC)=IQ(KQSP+LCDIR+LC1+KFLCYC)
         ENDIF
         IF(LCOLD.NE.0.AND.LCOLD.NE.LC1)THEN
            LC1=LCOLD
            GO TO 51
         ENDIF
      ENDIF
C
C           N option given. return neighbours
C
      IF(IOPTN.NE.0)THEN
         IF(I.EQ.1)THEN
            IQUEST(30)=0
         ELSE
            IQUEST(30)=NWKEY
            DO 52 J=1,NWKEY
               IF(J.LT.10)THEN
                  LKCJ=LK+(NWKEY+1)*(I-2)
                  IQUEST(30+J)=IQ(KQSP+LCDIR+LKCJ+J)
                  IKDES=(J-1)/10
                  IKBIT1=3*J-30*IKDES-2
                  IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).GE.3)THEN
                     CALL ZITOH(IQUEST(30+J),IQUEST(30+J),1)
                  ENDIF
               ENDIF
  52        CONTINUE
         ENDIF
         IF(I.EQ.NKEYS)THEN
            IQUEST(40)=0
         ELSE
            IQUEST(40)=NWKEY
            DO 53 J=1,NWKEY
               IF(J.LT.10)THEN
                  LKCJ=LK+(NWKEY+1)*I
                  IQUEST(40+J)=IQ(KQSP+LCDIR+LKCJ+J)
                  IKDES=(J-1)/10
                  IKBIT1=3*J-30*IKDES-2
                  IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).GE.3)THEN
                     CALL ZITOH(IQUEST(40+J),IQUEST(40+J),1)
                  ENDIF
               ENDIF
  53        CONTINUE
         ENDIF
      ENDIF
      GO TO 99
*
*           Error
*
  90  IQUEST(1)=1
      IF(IOPTN.NE.0)THEN
         IF(NKEYS.GT.0)THEN
            IQUEST(30)=NWKEY
            IQUEST(40)=NWKEY
            DO 91 J=1,NWKEY
               IF(J.GE.10)GO TO 91
               LKCJ=LK+(NWKEY+1)*(NKEYS-1)
               IQUEST(30+J)=IQ(KQSP+LCDIR+LK+J)
               IQUEST(40+J)=IQ(KQSP+LCDIR+LKCJ+J)
               IKDES=(J-1)/10
               IKBIT1=3*J-30*IKDES-2
               IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).GE.3)THEN
                  CALL ZITOH(IQUEST(30+J),IQUEST(30+J),1)
                  CALL ZITOH(IQUEST(40+J),IQUEST(40+J),1)
               ENDIF
  91        CONTINUE
         ENDIF
      ENDIF
*
  99  RETURN
      END
