*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:41  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZCHCY
*
************************************************************************
*
*      To check cycles block and reset pointers if needed
*
* Called by <USER>,RZINPA,RZLDIR
*
*  Author  :
*  Written :
*  Last mod: 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
*
*-----------------------------------------------------------------------
*
*
      IQUEST(1) = 0
      LD=IQ(KQSP+LCDIR+KLD)
      LK=IQ(KQSP+LCDIR+KLK)
      LF=IQ(KQSP+LCDIR+KLF)
      LC=IQ(KQSP+LCDIR+KLC)
      LE=IQ(KQSP+LCDIR+KLE)
      NKEYS =IQ(KQSP+LCDIR+KNKEYS)
      NWKEY =IQ(KQSP+LCDIR+KNWKEY)
      LC0   = LE-NKEYS*KLCYCL+1
      IF(LC.EQ.LC0) RETURN
*
*     Check consistancy of cycles block
*
      NGOOD = 0
      NBAD  = 0
      DO 10 LKC=LC,LE-KLCYCL+1,KLCYCL
         IF(IQ(KQSP+LCDIR+LKC).EQ.-1)THEN
            NBAD = NBAD + 1
         ELSE
            NGOOD = NGOOD + 1
            IF(NGOOD.EQ.1) LC1 = LKC
         ENDIF
   10 CONTINUE
      LC2 = NBAD*KLCYCL+LC
      IF(LC0.EQ.LC1.AND.LC0.EQ.LC2.AND.NGOOD.EQ.NKEYS) THEN
         IQ(KQSP+LCDIR+KLC)    = LC0
*
*     Mark directory as modified
*
         CALL SBIT1(IQ(KQSP+LCDIR),2)
         CALL SBIT1(IQ(KQSP+LTOP),2)
      ELSE
         IQUEST(1) = 1
      ENDIF
 
      END
