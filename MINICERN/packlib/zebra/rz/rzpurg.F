*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:05  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPURG(NKEEP)
*
************************************************************************
*
*        Purge cycles in the CWD
* Input:
*   NKEEP   Number of cycles which must be kept for the given key
*           If NKEEP < 1 then NKEEP is taken to be 1 and only the highest
*           cycle is kept
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 06.04.86
*  Last mod: 11.12.88
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      NK=NKEEP
      IF(NK.LT.1)NK=1
*
*           Check if write permission
*
      IF(LQRS.EQ.0)GO TO 999
      IFLAG=1
      CALL RZMODS('RZPURG',IFLAG)
      IF(IFLAG.NE.0)GO TO 999
*
      LK=IQ(KQSP+LCDIR+KLK)
      LC=IQ(KQSP+LCDIR+KLC)
      LE=IQ(KQSP+LCDIR+KLE)
      NKEYS =IQ(KQSP+LCDIR+KNKEYS)
      NWKEY =IQ(KQSP+LCDIR+KNWKEY)
      IF(NKEYS.EQ.0)GO TO 999
*
      NPUOLD=0
      IF(LPURG.NE.0)THEN
         NPURG=IQ(KQSP+LPURG+1)
         DO 5 I=1,NPURG
            NPUOLD=NPUOLD+IQ(KQSP+LPURG+2*I+1)-IQ(KQSP+LPURG+2*I)+1
   5     CONTINUE
      ENDIF
*
      DO 20 IK=1,NKEYS
         LKC=LK+(NWKEY+1)*(IK-1)
         LCYC=IQ(KQSP+LCDIR+LKC)
         NC=0
  10     NC=NC+1
         IF (KVSCYC.EQ.0) THEN
            LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC+KPPCYC), 1,16)
         ELSE
            LCOLD = IQ(KQSP+LCDIR+LCYC+KPPCYC)
         ENDIF
*
*              Check for first cycle to be kept
*
         LKEEP=LCOLD
         IF(LKEEP.NE.0)THEN
  12        IF(JBIT(IQ(KQSP+LCDIR+LKEEP+KFLCYC),5).EQ.0)THEN
               IF (KVSCYC.EQ.0) THEN
                  LKEEP = JBYT(IQ(KQSP+LCDIR+LKEEP+KPPCYC),1,16)
               ELSE
                  LKEEP = IQ(KQSP+LCDIR+LKEEP+KPPCYC)
               ENDIF
               IF(LKEEP.NE.0)GO TO 12
            ENDIF
         ENDIF
         IF(NC.EQ.NK)THEN
            IF (KVSCYC.EQ.0) THEN
               CALL SBYT(LKEEP,IQ(KQSP+LCDIR+LCYC+KPPCYC),1,16)
            ELSE
               IQ(KQSP+LCDIR+LCYC+KPPCYC) = LKEEP
            ENDIF
         ENDIF
         IF(NC.GT.NK)THEN
            IF(JBIT(IQ(KQSP+LCDIR+LCYC+KFLCYC),5).NE.0)THEN
               IF (KVSCYC.EQ.0) THEN
                  CALL SBYT(LKEEP,IQ(KQSP+LCDIR+LCYC+KPPCYC),1,16)
               ELSE
                  IQ(KQSP+LCDIR+LCYC+KPPCYC) = LKEEP
               ENDIF
               GO TO 15
            ENDIF
            IF (KVSCYC.EQ.0) THEN
               IR1  = JBYT(IQ(KQSP+LCDIR+LCYC+KFRCYC),17,16)
               IP1  = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,16)
               NW   = JBYT(IQ(KQSP+LCDIR+LCYC+KNWCYC), 1,20)
               IR2  = JBYT(IQ(KQSP+LCDIR+LCYC+KSRCYC),17,16)
            ELSE
               IR1  = IQ(KQSP+LCDIR+LCYC+KFRCYC)
               IP1  = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,20)
               NW   = IQ(KQSP+LCDIR+LCYC+KNWCYC)
               IR2  = IQ(KQSP+LCDIR+LCYC+KSRCYC)
            ENDIF
            IRL  =0
            NWL  =0
*
*              Mark all records that can be purged in first pass
*
            NLEFT=LREC-IP1+1
            NW1=NW
            IF(NW1.GE.NLEFT)NW1=NLEFT
            IF(IR2.NE.0)THEN
               NR=(NW-NW1-1)/LREC+1
               IF(NR.GT.1)THEN
                  CALL RZPURF(NR-1,IR2)
               ENDIF
               IRL=IR2+NR-1
               NWL=NW-NW1-(NR-1)*LREC
            ENDIF
            IF(NW1.EQ.LREC)THEN
               CALL RZPURF(1,IR1)
               IR1=0
            ENDIF
            IF(NWL.EQ.LREC)THEN
               CALL RZPURF(1,IRL)
               IRL=0
            ENDIF
            IRLOUT=IQ(KQSP+LCDIR+KRLOUT)
            IF(IRL.EQ.IRLOUT.AND.NWL+1.EQ.IQ(KQSP+LCDIR+KIP1))THEN
               CALL RZPURF(1,IRL)
            ENDIF
            IQ(KQSP+LCDIR+LCYC  )=-1
            IQ(KQSP+LCDIR+LCYC+1)=IR1
            IQ(KQSP+LCDIR+LCYC+2)=IRL
         ENDIF
*
  15     IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
            LCYC=LCOLD
            GO TO 10
         ENDIF
  20  CONTINUE
*
*           Now loop on all purged cycles to find complete records
*           purged
*
      DO 70 LKC=LC,LE-KLCYCL+1,KLCYCL
         IF(IQ(KQSP+LCDIR+LKC).NE.-1)GO TO 70
         IR1=IQ(KQSP+LCDIR+LKC+1)
         IRL=IQ(KQSP+LCDIR+LKC+2)
         IF(IR1.NE.0)THEN
            DO 30 LKC1=LC,LE-KLCYCL+1,KLCYCL
               IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)GO TO 30
               IF (KVSCYC.EQ.0) THEN
                  KR1  = JBYT(IQ(KQSP+LCDIR+LKC1+KFRCYC),17,16)
                  KP1  = JBYT(IQ(KQSP+LCDIR+LKC1+KORCYC), 1,16)
                  NW   = JBYT(IQ(KQSP+LCDIR+LKC1+KNWCYC), 1,20)
                  KR2  = JBYT(IQ(KQSP+LCDIR+LKC1+KSRCYC),17,16)
               ELSE
                  KR1  = IQ(KQSP+LCDIR+LKC1+KFRCYC)
                  KP1  = JBYT(IQ(KQSP+LCDIR+LKC1+KORCYC), 1,20)
                  NW   = IQ(KQSP+LCDIR+LKC1+KNWCYC)
                  KR2  = IQ(KQSP+LCDIR+LKC1+KSRCYC)
               ENDIF
               KRL  =0
               NLEFT=LREC-KP1+1
               NW1=NW
               IF(NW1.GE.NLEFT)NW1=NLEFT
               IF(KR2.NE.0)THEN
                  NR=(NW-NW1-1)/LREC+1
                  KRL=KR2+NR-1
               ENDIF
               IF(KR1.EQ.IR1.OR.KRL.EQ.IR1)GO TO 40
  30        CONTINUE
            CALL RZPURF(1,IR1)
         ENDIF
*
  40     IF(IRL.NE.0)THEN
            DO 50 LKC1=LC,LE-KLCYCL+1,KLCYCL
               IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)GO TO 50
               IF (KVSCYC.EQ.0) THEN
                  KR1  = JBYT(IQ(KQSP+LCDIR+LKC1+KFRCYC),17,16)
                  KP1  = JBYT(IQ(KQSP+LCDIR+LKC1+KORCYC), 1,16)
                  NW   = JBYT(IQ(KQSP+LCDIR+LKC1+KNWCYC), 1,20)
                  KR2  = JBYT(IQ(KQSP+LCDIR+LKC1+KSRCYC),17,16)
               ELSE
                  KR1  = IQ(KQSP+LCDIR+LKC1+KFRCYC)
                  KP1  = JBYT(IQ(KQSP+LCDIR+LKC1+KORCYC), 1,20)
                  NW   = IQ(KQSP+LCDIR+LKC1+KNWCYC)
                  KR2  = IQ(KQSP+LCDIR+LKC1+KSRCYC)
               ENDIF
               KRL  =0
               NLEFT=LREC-KP1+1
               NW1=NW
               IF(NW1.GE.NLEFT)NW1=NLEFT
               IF(KR2.NE.0)THEN
                  NR=(NW-NW1-1)/LREC+1
                  KRL=KR2+NR-1
               ENDIF
               IF(KR1.EQ.IRL.OR.KRL.EQ.IRL)GO TO 70
  50        CONTINUE
            CALL RZPURF(1,IRL)
         ENDIF
  70  CONTINUE
*
*           Garbage collection on cycles area + relocation
*
      LKC3=LE-KLCYCL+1
  80  IF(LKC3.LT.LC)GO TO 200
      IF(IQ(KQSP+LCDIR+LKC3).EQ.-1)THEN
         LKC3=LKC3+KLCYCL
         LKC2=LKC3-2*KLCYCL
  90     IF(IQ(KQSP+LCDIR+LKC2).NE.-1)THEN
            LKC2=LKC2+KLCYCL
            LKC1=LKC2-2*KLCYCL
            IF(LKC1.LT.LC)LKC1=LC
 100        IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)THEN
               IF(LKC1.GT.LC)LKC1=LKC1+KLCYCL
            ELSE
               IF(LKC1.GT.LC)THEN
                  LKC1=LKC1-KLCYCL
                  GO TO 100
               ENDIF
            ENDIF
*
            NPUSH=LKC3-LKC2
            DO 110 LKC=LC,LKC2-KLCYCL,KLCYCL
               IF(IQ(KQSP+LCDIR+LKC).NE.-1)THEN
                  IF (KVSCYC.EQ.0) THEN
                     LCOLD = JBYT(IQ(KQSP+LCDIR+LKC+KPPCYC),1,16)
                  ELSE
                     LCOLD = IQ(KQSP+LCDIR+LKC+KPPCYC)
                  ENDIF
                  IF(LCOLD.GE.LKC1.AND.LCOLD.LT.LKC2)THEN
                     LCOLD=LCOLD+NPUSH
                     IF (KVSCYC.EQ.0) THEN
                        CALL SBYT(LCOLD,IQ(KQSP+LCDIR+LKC+KPPCYC),1,16)
                     ELSE
                        IQ(KQSP+LCDIR+LKC+KPPCYC) = LCOLD
                     ENDIF
                  ENDIF
               ENDIF
 110        CONTINUE
*
            DO 120 IK=1,NKEYS
               LCYC=IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))
               IF(LCYC.GE.LKC1.AND.LCYC.LT.LKC2)THEN
                  IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))=
     +            IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))+NPUSH
               ENDIF
 120        CONTINUE
*
            CALL UCOPY2(IQ(KQSP+LCDIR+LKC1),IQ(KQSP+LCDIR+LKC1+NPUSH),
     +                  LKC2-LKC1)
            LKC3=LKC1+NPUSH
            LKC2=LKC1-KLCYCL
            IF(LKC1.NE.LC)GO TO 90
            LC=LC+NPUSH
            GO TO 200
*
         ELSE
            LKC2=LKC2-KLCYCL
            IF(LKC2.GE.LC)GO TO 90
         ENDIF
      ELSE
         LKC3=LKC3-KLCYCL
         GO TO 80
      ENDIF
*
*           Reset internal pointers
*
 200  CONTINUE
*
      NPUNEW=0
      IF(LPURG.NE.0)THEN
         NPURG=IQ(KQSP+LPURG+1)
         DO 210 I=1,NPURG
            NPUNEW=NPUNEW+IQ(KQSP+LPURG+2*I+1)-IQ(KQSP+LPURG+2*I)+1
 210     CONTINUE
      ENDIF
      NPU=NPUNEW-NPUOLD
      IQUEST( 9)=IQ(KQSP+LCDIR+KQUOTA)-IQ(KQSP+LCDIR+KRUSED)+NPU
      IQUEST(11)=(LC-IQ(KQSP+LCDIR+KLC))/KLCYCL
      IQUEST(12)=NPU*LREC
      IQUEST(13)=NPU
*
      IQ(KQSP+LCDIR+KRUSED)=IQ(KQSP+LCDIR+KRUSED)-NPU
      NWP=NPU*LREC
      NMEGA=NWP/1000000
      IQ(KQSP+LCDIR+KMEGA)=IQ(KQSP+LCDIR+KMEGA)-NMEGA
      NWP=NWP-1000000*NMEGA
      IQ(KQSP+LCDIR+KWUSED)=IQ(KQSP+LCDIR+KWUSED)-NWP
      IQ(KQSP+LCDIR+KNFREE)=IQ(KQSP+LCDIR+KNFREE)+LC-IQ(KQSP+LCDIR+KLC)
      IQ(KQSP+LCDIR+KLC)=LC
*
 999  RETURN
      END
