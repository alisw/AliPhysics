#ifndef CERNLIB_ISAJET_WWPAR_INC
#define CERNLIB_ISAJET_WWPAR_INC
*
*
* wwpar.inc
*
      COMMON/WWPAR/SWW,TWW,UWW,WM2,ZM2,P1WW(5),P2WW(5),P3WW(5),P4WW(5)
     $,PZERO(4,4),S13,P3(5),Q1(5),Q3(5),JQWW(2)
     $,CQ,CV,CA,CV1,CA1,CV3,CA3,CS,CT,CU
      SAVE /WWPAR/
#if defined(CERNLIB_DOUBLE)
#include "isajet/wwpar2.inc"
* Ignoring t=pass
#endif
#if defined(CERNLIB_SINGLE)
#include "isajet/wwpar1.inc"
* Ignoring t=pass
#endif
      INTEGER   JQWW

#endif
