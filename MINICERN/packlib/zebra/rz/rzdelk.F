*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:45  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZDELK(KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*       To delete one or all keys in the CWD
* Input:
*   KEYU    Key array of dimension NWKEY (see RZMDIR)
*   ICYCLE  Cycle number of the key to be deleted
*           ICYCLE > highest cycle number means delete the highest cycle
*           ICYCLE = 0 means delete the lowest cycle
*           ICYCLE = -1, -2,... means delete the highest cycle -1, -2,...
*   CHOPT   Character variable specifying the options selected.
*           default
*                 Delete the explicitly specified cycle ICYCLE only.
*                 If cycle ICYCLE does not exist, no action is taken.
*           'C'   Delete  ALL cycles corresponding to  key (ICYCLE not used)
*           'S'   Delete all cycles smaller then cycle ICYCLE.
*           'K'   Delete  ALL cycles for all Keys (KEYU,ICYCLE not used )
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 20.04.86
*  Last mod: 16.04.93 JDS. Return codes, deletion of objects at beginning
*                          of cycles structure
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*  IQUEST(1) = 0: ok
*  IQUEST(1) = 1: specified object not found
*  IQUEST(1) = 2: directory is empty
*  IQUEST(1) = 3: no RZ control bank (LQRS=0)
*  IQUEST(1) = 4: no write permission
*  IQUEST(1) =11: key/cycle discrepency
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*(*) CHOPT
      DIMENSION KEYU(*)
      DIMENSION    IOPTV(3)
      EQUIVALENCE (IOPTC,IOPTV(1)), (IOPTS,IOPTV(2))
      EQUIVALENCE (IOPTK,IOPTV(3))
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
      IQUEST(1) = 0
      IQWARN    = 0
      CALL UOPTC(CHOPT,'CSK',IOPTV)
*
*           Check if write permission
*
      IF(LQRS.EQ.0) THEN
         IQUEST(1) = 3
         GOTO 230
      ENDIF
 
      IFLAG=1
      CALL RZMODS('RZDELK',IFLAG)
      IF(IFLAG.NE.0) THEN
         IQUEST(1) = 4
         GOTO 230
      ENDIF
*
      LD=IQ(KQSP+LCDIR+KLD)
      LK=IQ(KQSP+LCDIR+KLK)
      LF=IQ(KQSP+LCDIR+KLF)
      LC=IQ(KQSP+LCDIR+KLC)
      LE=IQ(KQSP+LCDIR+KLE)
      NKEYS =IQ(KQSP+LCDIR+KNKEYS)
      NWKEY =IQ(KQSP+LCDIR+KNWKEY)
      IF(NKEYS.EQ.0)GOTO 220
*
*     Look for cycles marked for deletion by a previous call
*     Set IQUEST(2) to warn application to rebuild bit map
*
      DO 10 LKC=LC,LE-KLCYCL+1,KLCYCL
         IF(IQ(KQSP+LCDIR+LKC).EQ.-1)GOTO 20
   10 CONTINUE
      GOTO 30
   20 IQWARN = 1
   30 CONTINUE
*
      NPUOLD=0
      IF(LPURG.NE.0)THEN
         NPURG=IQ(KQSP+LPURG+1)
         DO 40 I=1,NPURG
            NPUOLD=NPUOLD+IQ(KQSP+LPURG+2*I+1)-IQ(KQSP+LPURG+2*I)+1
   40    CONTINUE
      ENDIF
*
*           Option K. delete all keys
*
      IF(IOPTK.NE.0)THEN
         DO 50 LKC=LC,LE-KLCYCL+1,KLCYCL
            IF (KVSCYC.EQ.0) THEN
               IR1 = JBYT(IQ(KQSP+LCDIR+LKC+KFRCYC),17,16)
               IP1 = JBYT(IQ(KQSP+LCDIR+LKC+KORCYC), 1,16)
               NW  = JBYT(IQ(KQSP+LCDIR+LKC+KNWCYC), 1,20)
               IR2 = JBYT(IQ(KQSP+LCDIR+LKC+KSRCYC),17,16)
            ELSE
               IR1 = IQ(KQSP+LCDIR+LKC+KFRCYC)
               IP1 = JBYT(IQ(KQSP+LCDIR+LKC+KORCYC), 1,20)
               NW  = IQ(KQSP+LCDIR+LKC+KNWCYC)
               IR2 = IQ(KQSP+LCDIR+LKC+KSRCYC)
            ENDIF
            NLEFT=LREC-IP1+1
            IF(NW.LE.NLEFT)THEN
               NR=0
            ELSE
               NR=(NW-NLEFT-1)/LREC+1
            ENDIF
            IF(IR2.EQ.IR1+1)THEN
               CALL RZPURF(NR+1,IR1)
            ELSE
               CALL RZPURF(1,IR1)
               IF(NR.NE.0)CALL RZPURF(NR,IR2)
            ENDIF
   50    CONTINUE
         LF=LK
         LC=LE+1
         NRD=IQ(KQSP+LCDIR+LD)
         IQ(KQSP+LCDIR+KRUSED) =NRD
         IQ(KQSP+LCDIR+KWUSED) =NRD*LREC
         IQ(KQSP+LCDIR+KMEGA)  =0
         IQ(KQSP+LCDIR+KIP1)   =1
         IQ(KQSP+LCDIR+KNFREE) =LC-LF
         IQ(KQSP+LCDIR+KLF)    =LF
         IQ(KQSP+LCDIR+KLC)    =LC
         IQ(KQSP+LCDIR+KNKEYS) =0
         GOTO 240
      ENDIF
*
*           Search KEY and CYCLE
*
      DO 60 I=1,NWKEY
         IKDES=(I-1)/10
         IKBIT1=3*I-30*IKDES-2
         IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
            KEY(I)=KEYU(I)
         ELSE
            CALL ZHTOI(KEYU(I),KEY(I),1)
         ENDIF
   60 CONTINUE
      DO 80 I=1,NKEYS
         DO 70 K=1,NWKEY
            LKC=LK+(NWKEY+1)*(I-1)
            IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GOTO 80
   70    CONTINUE
         LCYC  =IQ(KQSP+LCDIR+LKC)
         LCPRE =LCYC
         LKK   =LKC
         IF (KVSCYC.NE.0) THEN
*           IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.I) GO TO 250
            IF (IQ(KQSP+LCDIR+LCYC+KKYCYC).NE.
     +          IQ(KQSP+LCDIR+LKC+1)) GO TO 250
         ENDIF
         ICTOP =JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
         KCYCLE=ICYCLE
         IF(KCYCLE.GT.ICTOP)KCYCLE=ICTOP
         GOTO 90
   80 CONTINUE
      GOTO 210
*
*              Do we keep this cycle ?
*
   90 IF (KVSCYC.EQ.0) THEN
         LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC+KPPCYC), 1,16)
      ELSE
         LCOLD = IQ(KQSP+LCDIR+LCYC+KPPCYC)
      ENDIF
      IDEL=0
      ICY =JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
      IF(KCYCLE.EQ.ICY.AND.IOPTS.EQ.0)IDEL=1
      IF(IOPTC.NE.0)IDEL=1
      IF(IOPTK.NE.0)IDEL=1
      IF(IOPTS.NE.0)THEN
         IF(ICY.LT.ICYCLE)IDEL=1
      ENDIF
      IF(ICYCLE.EQ.0.AND.LCOLD.EQ.0)IDEL=1
      IF(ICYCLE.LT.0)THEN
         IF(ICY.EQ.ICTOP-ICYCLE)IDEL=1
      ENDIF
*
*              Mark all records that can be purged in first pass
*
      IF(IDEL.NE.0)THEN
         IF(ICY.EQ.ICTOP)THEN
            IQ(KQSP+LCDIR+LKK)=LCOLD
         ELSE
            IF(LCOLD.EQ.0.AND.IOPTC.NE.0)THEN
               IQ(KQSP+LCDIR+LKK)=0
               IQ(KQSP+LCDIR+LCPRE)=-1
            ELSE
               CALL SBYT(LCOLD,IQ(KQSP+LCDIR+LCPRE),1,16)
            ENDIF
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
         NLEFT=LREC-IP1+1
         NW1=NW
         IF(NW1.GE.NLEFT)NW1=NLEFT
         IF(IR2.NE.0)THEN
            NR=(NW-NW1-1)/LREC+1
            IF(NR.GT.1) CALL RZPURF(NR-1,IR2)
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
            IRL=0
         ENDIF
         IQ(KQSP+LCDIR+LCYC)=-1
         IQ(KQSP+LCDIR+LCYC+1)=IR1
         IQ(KQSP+LCDIR+LCYC+2)=IRL
         IQ(KQSP+LCDIR+LCYC+3)=NWL
         IF(IR1.NE.0)CALL SBYT(NW1,IQ(KQSP+LCDIR+LCYC+3),21,12)
      ELSE
         LCPRE=LCYC
      ENDIF
*
      IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
         IF(KCYCLE.LT.ICY.OR.IOPTS.NE.0.OR.IOPTK.NE.0.OR.IOPTC.NE.0)THEN
            LCYC=LCOLD
            GOTO 90
         ENDIF
      ENDIF
*
*           Now loop on all purged cycles to find complete records
*           purged
*
      DO 130 LKC=LC,LE-KLCYCL+1,KLCYCL
         IF(IQ(KQSP+LCDIR+LKC).NE.-1)GOTO 130
         IR1=IQ(KQSP+LCDIR+LKC+1)
         IRL=IQ(KQSP+LCDIR+LKC+2)
         IF(IR1.NE.0)THEN
            DO 100 LKC1=LC,LE-KLCYCL+1,KLCYCL
               IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)GOTO 100
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
               IF(KR1.EQ.IR1.OR.KRL.EQ.IR1)GOTO 110
  100       CONTINUE
            CALL RZPURF(1,IR1)
         ENDIF
*
  110    IF(IRL.NE.0)THEN
            DO 120 LKC1=LC,LE-KLCYCL+1,KLCYCL
               IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)GOTO 120
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
               IF(KR1.EQ.IRL.OR.KRL.EQ.IRL)GOTO 130
  120       CONTINUE
            CALL RZPURF(1,IRL)
         ENDIF
  130 CONTINUE
*
*           Garbage collection on cycles area + relocation
*
      LKC3=LE-KLCYCL+1
  140 IF(LKC3.LT.LC)GOTO 190
*
*     Found a deleted object. Now look for previous undeleted object
*
      IF(IQ(KQSP+LCDIR+LKC3).EQ.-1)THEN
 
         LKC3=LKC3+KLCYCL
         LKC2=LKC3-2*KLCYCL
*
*     First object?
*
         IF(LKC2.LT.LC) THEN
            LC    = LKC3
            GOTO 190
         ENDIF
 
  150    IF(IQ(KQSP+LCDIR+LKC2).NE.-1)THEN
            LKC2=LKC2+KLCYCL
            LKC1=LKC2-2*KLCYCL
            IF(LKC1.LT.LC)LKC1=LC
  160       IF(IQ(KQSP+LCDIR+LKC1).EQ.-1)THEN
                   IF(LKC1.GT.LC .OR.
     +             (LKC1.EQ.LC.AND.LKC2.NE.LC))LKC1=LKC1+KLCYCL
            ELSE
               IF(LKC1.GT.LC)THEN
                  LKC1=LKC1-KLCYCL
                  GOTO 160
               ENDIF
            ENDIF
*
            NPUSH=LKC3-LKC2
*
*      Update pointers in cycles block
*
            DO 170 LKC=LC,LKC2-KLCYCL,KLCYCL
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
  170       CONTINUE
*
*      Update pointers from KEYS block to CYCLES block
*
            DO 180 IK=1,NKEYS
               LCYC=IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))
               IF(LCYC.GE.LKC1.AND.LCYC.LT.LKC2)THEN
                  IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))=
     +            IQ(KQSP+LCDIR+LK+(NWKEY+1)*(IK-1))+NPUSH
               ENDIF
  180       CONTINUE
*
*     Squeeze out deleted cycles
*
            CALL UCOPY2(IQ(KQSP+LCDIR+LKC1),IQ(KQSP+LCDIR+LKC1+NPUSH),
     +                  LKC2-LKC1)
            LKC3=LKC1+NPUSH
            LKC2=LKC1-KLCYCL
*
*     Only deleted objects before this block?
*
            IF(IQ(KQSP+LCDIR+LKC2).EQ.-1.AND.LKC2.EQ.LC) THEN
               LC  = LKC3
               GOTO 190
            ENDIF
 
            IF(LKC1.NE.LC)GOTO 150
            LC=LC+NPUSH
            GOTO 190
*
         ELSE
            LKC2=LKC2-KLCYCL
*
*     Only deleted objects before this block?
*
            IF(IQ(KQSP+LCDIR+LKC2).EQ.-1.AND.LKC2.EQ.LC) THEN
               LC  = LKC3
               GOTO 190
            ENDIF
 
            IF(LKC2.GE.LC)GOTO 150
 
         ENDIF
      ELSE
         LKC3=LKC3-KLCYCL
         GOTO 140
      ENDIF
*
*           Remove KEY from K area if only one cycle
*
  190 CONTINUE
      IF(IQ(KQSP+LCDIR+LKK).EQ.0)THEN
         IF(LKK+NWKEY+1.LT.LF)THEN
            CALL UCOPY2(IQ(KQSP+LCDIR+LKK+NWKEY+1),
     +                  IQ(KQSP+LCDIR+LKK),LF-LKK)
         ENDIF
         LF=LF-NWKEY-1
         IQ(KQSP+LCDIR+KNFREE)=IQ(KQSP+LCDIR+KNFREE)+NWKEY+1
         IQ(KQSP+LCDIR+KLF)=LF
         IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)-1
         IF(IQ(KQSP+LCDIR+KNKEYS).EQ.0)LC=LE+1
      ENDIF
*
*           Reset internal pointers
*
      NPUNEW=0
      IF(LPURG.NE.0)THEN
         NPURG=IQ(KQSP+LPURG+1)
         DO 200 I=1,NPURG
            NPUNEW=NPUNEW+IQ(KQSP+LPURG+2*I+1)-IQ(KQSP+LPURG+2*I)+1
  200    CONTINUE
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
      GOTO 240
*
*           KEY not found
*
  210 IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,10000)
10000 FORMAT(' RZDELK. Key not found')
      IQUEST(1)=1
      IQUEST(2)=IQWARN
      RETURN
 
  220 IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,10100)
10100 FORMAT(' RZDELK. directory is empty')
      IQUEST(1)=2
      IQUEST(2)=IQWARN
      RETURN
*
  230 IQUEST(2)=IQWARN
      RETURN
 
  240 IQUEST(2)=IQWARN
      RETURN
 
  250 IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,10100)
10200 FORMAT(' RZDELK. mismatch in key/cycle pointing')
      IQUEST(1)=11
      IQUEST(2)=IQWARN
      RETURN
      END
