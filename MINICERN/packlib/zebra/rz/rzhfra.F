*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:52  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
       SUBROUTINE RZHFRA(IH,NW)
*
***********************************************************************
*
*       Converts array of NW ASCII words into  EBCDIC hollerith
*
*  Author  : R.Brun CN/AS
*  Written : 17/10/92
*  Last mod: 17/10/92
*
***********************************************************************
#if defined(CERNLIB_QMIBM)
      DIMENSION IH(1)
      PARAMETER          (NCHAR = 126)
      DIMENSION IEB(NCHAR)
      DATA IEB /8*0,5,22*0,
     + 64, 90,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, 75, 97,
     +240,241,242,243,244,245,246,247,248,249,122, 94, 76,126,110,111,
     +124,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214,
     +215,216,217,226,227,228,229,230,231,232,233,173,224,189, 95,109,
     +121,129,130,131,132,133,134,135,136,137,145,146,147,148,149,150,
     +151,152,153,162,163,164,165,166,167,168,169,192, 79,208,161/
*-----------------------------------
*

#include "zebra/q_jbyt.inc"

      DO 20 IW=1,NW
         DO 10 I=1,25,8
            IBYTE=JBYT(IH(IW),I,8)
            IF(IBYTE.GT.126)IBYTE=32
            CALL SBYT(IEB(IBYTE),IH(IW),I,8)
  10     CONTINUE
  20  CONTINUE
#endif
      END
