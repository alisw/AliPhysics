#ifndef CERNLIB_ISAJET_GETJET_INC
#define CERNLIB_ISAJET_GETJET_INC
*
*
* getjet.inc
*
      INTEGER   NJMAX
      PARAMETER (NJMAX=50)
      COMMON/GETCOM/JETNO(NCY,NCPHI),NCJET,PCJET(4,NJMAX),ETJET(NJMAX)
      SAVE /GETCOM/
#if defined(CERNLIB_LEVEL2)
#include "isajet/l2getj.inc"
* Ignoring t=pass
#endif
      INTEGER   JETNO,NCJET
      REAL      PCJET,ETJET

#endif
