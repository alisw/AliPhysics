#ifndef CERNLIB_ISAJET_DKYTAB_INC
#define CERNLIB_ISAJET_DKYTAB_INC
*
*
* dkytab.inc
*
C          LOOK must be dimensioned to the maximum value of INDEX.
      INTEGER   MXLOOK
      PARAMETER (MXLOOK=500)
      INTEGER   MXDKY
      PARAMETER (MXDKY=3000)
      COMMON/DKYTAB/LOOK(MXLOOK),CBR(MXDKY),MODE(5,MXDKY),MELEM(MXDKY)
      SAVE /DKYTAB/
#if defined(CERNLIB_LEVEL2)
#include "isajet/l2dky.inc"
* Ignoring t=pass
#endif
      INTEGER   LOOK,MODE,MELEM
      REAL      CBR

#endif
