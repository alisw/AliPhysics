*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:53  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
       SUBROUTINE RZHTOA(IH,NW)
*
***********************************************************************
*
*       Converts array of NW EBCDIC hollerith words into ASCII
*
*  Author  : R.Brun CN/AS
*  Written : 17/10/92
*  Last mod: 17/10/92
*
***********************************************************************
#if defined(CERNLIB_QMIBM)
      DIMENSION IH(1)
      PARAMETER          (NCHAR = 256)
      DIMENSION IASCI(NCHAR)
      DATA IASCI/
     +4*32,9,58*32,32,10*32,46,60,40,43,124,38,9*32,33,36,42,41,59,94,
     +45,47,9*32,44,37,95,62,63,9*32,96,58,35,64,39,61,34,32,97,98,99,
     +100,101,102,103,104,105,7*32,106,107,108,109,110,111,112,113,114,
     +7*32,126,115,116,117,118,119,120,121,122,3*32,91,15*32,93,32,32,
     +123,65,66,67,68,69,70,71,72,73,6*32,125,74,75,76,77,78,79,80,81,
     +82,6*32,92,32,83,84,85,86,87,88,89,90,6*32,48,49,50,51,52,53,54,
     +55,56,57,7*32/
*-----------------------------------
*

#include "zebra/q_jbyt.inc"

      DO 20 IW=1,NW
         DO 10 I=1,25,8
            IBYTE=JBYT(IH(IW),I,8)
            CALL SBYT(IASCI(IBYTE),IH(IW),I,8)
  10     CONTINUE
  20  CONTINUE
#endif
      END
