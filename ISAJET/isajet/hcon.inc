#ifndef CERNLIB_ISAJET_HCON_INC
#define CERNLIB_ISAJET_HCON_INC
*
*
* hcon.inc
*
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
#if defined(CERNLIB_DOUBLE)
#include "isajet/hcon2.inc"
* Ignoring t=pass
#endif
#if defined(CERNLIB_SINGLE)
#include "isajet/hcon1.inc"
* Ignoring t=pass
#endif
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS

#endif
