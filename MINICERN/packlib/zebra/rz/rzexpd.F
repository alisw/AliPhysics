*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:48  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZEXPD(CHROUT,NMORE)
*
************************************************************************
*
*          Routine to expand one directory by NMORE words
*
* Called by RZCOP1,RZMDIR,RZOUT,RZQUOT,RZVOUT
*
*  Author  : R.Brun DD/US/PD
*  Written : 05.04.86
*  Last mod: 18.06.92 - bug fix for the case when >1 record is allocated
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*6 CHROUT
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
      LD = IQ(KQSP+LCDIR+KLD)
      LB = IQ(KQSP+LCDIR+KLB)
      LS = IQ(KQSP+LCDIR+KLS)
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
      LC = IQ(KQSP+LCDIR+KLC)
      LE = IQ(KQSP+LCDIR+KLE)
      NPR= (NMORE-1)/LREC +1
      NPUSH=NPR*LREC
      CALL RZALLO(CHROUT,NPR,IALLOC)
      IF(IALLOC.EQ.0) GO TO 99
      CALL MZPUSH(JQPDVS,LCDIR,0,NPUSH,' ')
      NWFREE=IQ(KQSP+LCDIR+KNFREE)+NPUSH-NPR
      IQ(KQSP+LCDIR+KNFREE)=NWFREE
*
*           Move cycles
*
      IF(LC.LT.LE)THEN
         NKEYS=IQ(KQSP+LCDIR+KNKEYS)
         NWKEY=IQ(KQSP+LCDIR+KNWKEY)
         DO 20 I=1,NKEYS
            LKC=LK+(NWKEY+1)*(I-1)
            IQ(KQSP+LCDIR+LKC)=IQ(KQSP+LCDIR+LKC)+NPUSH
  20     CONTINUE
         DO 30 LKC=LC,LE-KLCYCL+1,KLCYCL
            IF (KVSCYC.EQ.0) THEN
               LCOLD = JBYT(IQ(KQSP+LCDIR+LKC+KPPCYC),1,16)
            ELSE
               LCOLD = IQ(KQSP+LCDIR+LKC+KPPCYC)
            ENDIF
            IF(LCOLD.NE.0)THEN
               LCOLD=LCOLD+NPUSH
               IF (KVSCYC.EQ.0) THEN
                  CALL SBYT(LCOLD,IQ(KQSP+LCDIR+LKC+KPPCYC),1,16)
               ELSE
                  IQ(KQSP+LCDIR+LKC+KPPCYC) = LCOLD
               ENDIF
            ENDIF
  30     CONTINUE
         CALL UCOPY2(IQ(KQSP+LCDIR+LC),
     +               IQ(KQSP+LCDIR+LC+NPUSH),LE-LC+1)
      ENDIF
      LC=LC+NPUSH
      LE=LE+NPUSH
*
*           Insert new record in the list of records for the CWD
*           Move B S and K blocks
*
      CALL UCOPY2(IQ(KQSP+LCDIR+LB),IQ(KQSP+LCDIR+LB+NPR),LF-LB)
      LB=LB+NPR
      LS=LS+NPR
      LK=LK+NPR
      LF=LF+NPR
      NRD=IQ(KQSP+LCDIR+LD)
      IQ(KQSP+LCDIR+LD)=NRD+NPR
      DO 40 I=1,NPR
         IQ(KQSP+LCDIR+LD+NRD+I)=IALLOC+I-1
  40  CONTINUE
      IQ(KQSP+LCDIR+KRUSED)=IQ(KQSP+LCDIR+KRUSED)+NPR
      IQ(KQSP+LCDIR+KWUSED)=IQ(KQSP+LCDIR+KWUSED)+NPUSH
      CALL RZUSED(NPR,IALLOC)
      IQ(KQSP+LCDIR+KLB)=LB
      IQ(KQSP+LCDIR+KLS)=LS
      IQ(KQSP+LCDIR+KLK)=LK
      IQ(KQSP+LCDIR+KLF)=LF
      IQ(KQSP+LCDIR+KLC)=LC
      IQ(KQSP+LCDIR+KLE)=LE
*
  99  RETURN
      END
