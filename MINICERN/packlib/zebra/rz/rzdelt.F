*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:46  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZDELT(CHDIR)
*
************************************************************************
*
*       To delete the tree CHDIR in the CWD
* Input:
*   CHDIR   Character  variable  specifying  the directory  name  of  the
*           subtree of the CWD.
*           ' ' means delete the complete CWD tree
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 22.04.86
*  LAST MOD: 09.01.91
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*(*) CHDIR
      DIMENSION IHDIR(4),ISD(NLPATM),IRD(NLPATM),NSD(NLPATM)
      LOGICAL RZSAME
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
      NCD=LEN(CHDIR)
      IF(NCD.GT.16)NCD=16
      CALL VBLANK(IHDIR,4)
      CALL UCTOH(CHDIR,IHDIR,4,NCD)
      CALL ZHTOI(IHDIR,IHDIR,4)
*
      LS = IQ(KQSP+LCDIR+KLS)
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
*
*           Check permission
*
      IFLAG=1
      CALL RZMODS('RZDELT',IFLAG)
      IF(IFLAG.NE.0)GO TO 99
*
*                Check if subdirectory exists
*
      NSDIR=IQ(KQSP+LCDIR+KNSD)
      IF(NSDIR.GT.0)THEN
         LRZ=LQ(KQSP+LCDIR-1)
         DO 20 I=1,NSDIR
            IF(RZSAME(IHDIR,IQ(KQSP+LCDIR+LS+7*(I-1)),4))THEN
               IOLD=LS+7*(I-1)
               IF (KVSCYC.EQ.0) THEN
                  IR1 = JBYT(IQ(KQSP+LCDIR+IOLD+5),1,18)
               ELSE
                  IR1 = IQ(KQSP+LCDIR+IOLD+5)
               ENDIF
               GO TO 25
            ENDIF
  20     CONTINUE
      ENDIF
      IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,1000)CHDIR
 1000 FORMAT(' RZDELT. Non existing directory, ',A)
      IQUEST(1)=1
      GO TO 99
*
*            If directory to be deleted is in memory, then
*            delete the corresponding tree
*
  25  IF(LRZ.NE.0)THEN
         IF(.NOT.RZSAME(IHDIR,IQ(KQSP+LRZ+1),4))THEN
            LRZ=LQ(KQSP+LRZ)
            GO TO 25
         ELSE
            CALL MZDROP(JQPDVS,LRZ,' ')
         ENDIF
      ENDIF
*
*           Remove directory name from D
*           Move K
*
      NWFREE=IQ(KQSP+LCDIR+KNFREE)+7
      CALL UCOPY2(IQ(KQSP+LCDIR+IOLD+7),IQ(KQSP+LCDIR+IOLD),LF-IOLD-7)
      LK=LK-7
      LF=LF-7
      NSDIR=NSDIR-1
      IQ(KQSP+LCDIR+KNFREE)=NWFREE
      IQ(KQSP+LCDIR+KNSD)=NSDIR
      IQ(KQSP+LCDIR+KLK)=LK
      IQ(KQSP+LCDIR+KLF)=LF
*
*            LRIN will be used as delete buffer.
*            Make sure it exists.
*
      IF(LRIN.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LRIN,LTOP,-7,'RZIN',0,0,LREC+1,2,-1)
      ENDIF
      IQ(KQSP+LTOP+KIRIN)=0
*
      CALL VZERO(ISD(2),19)
      NLEVEL=1
      ISD(1)=1
      IRD(1)=IR1
      NSD(1)=1
*
*            Read directory into buffer
*
  30  CALL RZIODO(LUN,LREC,IRD(NLEVEL),IQ(KQSP+LRIN+1),1)
      IF(IQUEST(1).NE.0)GO TO 90
      NSD(NLEVEL+1)=IQ(KQSP+LRIN+23)
      LDC=IQ(KQSP+LRIN+KLD)
      LCC=IQ(KQSP+LRIN+KLC)
      LEC=IQ(KQSP+LRIN+KLE)
      NRD=IQ(KQSP+LRIN+LDC)
      NPUSH=NRD*LREC-IQ(KQSP+LRIN-1)
      IF(NPUSH.GT.0)CALL MZPUSH(JQPDVS,LRIN,0,NPUSH,'I')
      IF(NRD.GT.1)THEN
         DO 40 I=1,NRD
            L1=KQSP+LRIN+(I-1)*LREC+1
            CALL RZIODO(LUN,LREC,IQ(KQSP+LRIN+LDC+I),IQ(L1),1)
            IF(IQUEST(1).NE.0)GO TO 90
  40     CONTINUE
      ENDIF
*
*            Delete all KEYS for this directory
*
      DO 50 LKC=LCC,LEC-KLCYCL+1,KLCYCL
         IF (KVSCYC.EQ.0) THEN
            IR1 = JBYT(IQ(KQSP+LRIN+LKC+KFRCYC),17,16)
            IP1 = JBYT(IQ(KQSP+LRIN+LKC+KORCYC), 1,16)
            NW  = JBYT(IQ(KQSP+LRIN+LKC+KNWCYC), 1,20)
            IR2 = JBYT(IQ(KQSP+LRIN+LKC+KSRCYC),17,16)
         ELSE
            IR1 = IQ(KQSP+LRIN+LKC+KFRCYC)
            IP1 = JBYT(IQ(KQSP+LRIN+LKC+KORCYC), 1,20)
            NW  = IQ(KQSP+LRIN+LKC+KNWCYC)
            IR2 = IQ(KQSP+LRIN+LKC+KSRCYC)
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
  50  CONTINUE
      DO 60 I=1,NRD
         CALL RZPURF(1,IQ(KQSP+LRIN+LDC+I))
  60  CONTINUE
*
*            Now look levels down
*
      NLEVEL=NLEVEL+1
  70  ISD(NLEVEL)=ISD(NLEVEL)+1
      IF(ISD(NLEVEL).LE.NSD(NLEVEL))THEN
         IS=ISD(NLEVEL)
         LSC=IQ(KQSP+LRIN+KLS)
         IF (KVSCYC.EQ.0) THEN
            IRD(NLEVEL) = JBYT(IQ(KQSP+LRIN+LSC+7*(IS-1)+5),1,18)
         ELSE
            IRD(NLEVEL) = IQ(KQSP+LRIN+LSC+7*(IS-1)+5)
         ENDIF
         GO TO 30
      ELSE
         ISD(NLEVEL)=0
         NLEVEL=NLEVEL-1
         IF(NLEVEL.GT.1)THEN
            CALL RZIODO(LUN,LREC,IRD(NLEVEL-1),IQ(KQSP+LRIN+1),1)
            IF(IQUEST(1).NE.0)GO TO 90
            GO TO 70
         ENDIF
      ENDIF
*
  90  CALL MZDROP(JQPDVS,LRIN,' ')
      LRIN=0
#if defined(CERNLIB_QMVAX)
         IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
*
  99  RETURN
      END
