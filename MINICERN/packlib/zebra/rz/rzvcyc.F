*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVCYC(LTAD)
*
************************************************************************
*
*           Fills up the offset parameters of cycle information
* Input:
*   LTAD    Address of the top directory
*
* Called by RZFDIR,RZFILE,RZMAKE
*
*  Author  : S.Banerjee
*  Written : 03.03.94
*  Last mod:
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
*-----------------------------------------------------------------------
      IF (LTAD.EQ.0) GO TO 99
      IF (IQ(KQSP+LTAD+KRZVER).EQ.0) THEN
         KLCYCL = 4
         KPPCYC = 0
         KFRCYC = 2
         KSRCYC = 0
         KFLCYC = 1
         KORCYC = 2
         KCNCYC = 3
         KNWCYC = 3
         KKYCYC =-1
         KVSCYC = 0
      ELSE
         KLCYCL = 7
         KPPCYC = 0
         KFRCYC = 2
         KSRCYC = 5
         KFLCYC = 1
         KORCYC = 3
         KCNCYC = 3
         KNWCYC = 4
         KKYCYC = 6
         KVSCYC = 1
      ENDIF
  99  RETURN
      END
